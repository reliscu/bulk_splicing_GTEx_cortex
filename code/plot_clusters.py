"""
plot_clusters.py

Genome-browser-style visualization of one splicing cluster.
Left panel  : exon/intron structure per event.
Right panel : per-cell-type r-value heatmap (from signif_info).

Usage
-----
from plot_clusters import plot_cluster

plot_cluster(cid, members[cid], exon_info, signif_info)
plot_cluster(cid, members[cid], exon_info, signif_info, save_path='cluster_5178.pdf')

signif_info format
------------------
signif_info[event_id]['gene_name'] = str                                  # store once per event
signif_info[event_id][cell_type]   = {'r': float, 'fdr': float, 'is_specific': bool}

Build it like this:
    for idx, row in signif_exons_df.iterrows():
        signif_info[idx]['gene_name'] = row['Gene']
        signif_info[idx][ct] = {'is_specific': row['is_specific'],
                                 'r': row['r'], 'fdr': row['fdr']}
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm


# ── shared row builder (both panels use the same sorted order) ─────────────
def _build_rows(event_ids, exon_info, signif_set):
    rows = []
    for ev in event_ids:
        if ev not in exon_info:
            continue
        meta = exon_info[ev]['meta']
        rows.append({
            'event': ev,
            'chrom': meta['chrom'],
            'gene':  meta.get('gene', ''),
            'es':    meta['es'],
            'ee':    meta['ee'],
            'us':    meta['us_intron_start'],
            'ds':    meta['ds_intron_end'],
            'signif': ev in signif_set,
        })
    rows.sort(key=lambda r: (r['es'], r['ee'], r['ds']))
    return rows


# ── left panel: exon / intron structure ───────────────────────────────────
def _plot_structure(ax, cid, rows, signif_set, signif_info):
    if not rows:
        ax.set_visible(False)
        return

    x_min = min(r['us'] for r in rows)
    x_max = max(r['ds'] for r in rows)
    span  = max(x_max - x_min, 1)
    pad   = span * 0.04

    for i, r in enumerate(rows):
        color = 'tomato' if r['signif'] else 'steelblue'

        # intron span line
        ax.plot([r['us'], r['ds']], [i, i],
                color='#bbbbbb', lw=0.8, zorder=1, solid_capstyle='round')

        # cassette exon block
        width = max(r['ee'] - r['es'], span * 0.003)
        ax.add_patch(Rectangle((r['es'], i - 0.3), width, 0.6,
                                facecolor=color, edgecolor='none', zorder=2))

        # event label
        parts  = r['event'].split('_')
        label  = '_'.join(parts[-2:]) if len(parts) >= 2 else r['event']
        marker = ' ★' if r['signif'] else ''
        ax.text(x_min - pad, i, label + marker,
                ha='right', va='center', fontsize=5.5,
                color='tomato' if r['signif'] else '#555555')

    # gene display name: scan all rows for first event that has gene_name in signif_info
    gene_display = next(
        (signif_info[r['event']].get('gene_name')
        for r in rows
        if r['event'] in signif_info and signif_info[r['event']].get('gene_name')),
        rows[0]['gene']
    )
    n_sig = sum(r['signif'] for r in rows)
    ax.set_title(
        f"cluster {cid}  |  {gene_display}  |  {rows[0]['chrom']}\n"
        f"{len(rows)} events  ({n_sig} significant in ≥1 cell type)",
        fontsize=6.5, loc='left', pad=3
    )
    ax.set_xlim(x_min - pad * 2, x_max + pad)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_yticks([])

    def _kb(x, _):
        return f'+{(x - x_min) / 1000:.1f}kb' if x != x_min else '0'
    ax.xaxis.set_major_formatter(plt.FuncFormatter(_kb))
    ax.xaxis.set_major_locator(plt.MaxNLocator(4, integer=False))
    ax.tick_params(axis='x', labelsize=5.5, pad=1)
    for spine in ('top', 'left', 'right'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)


# ── right panel: per-cell-type r-value heatmap ────────────────────────────
def _plot_celltype_panel(ax, rows, signif_info, cell_types):
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    cmap = plt.get_cmap('RdBu_r')

    for i, r in enumerate(rows):
        ev = r['event']
        for j, ct in enumerate(cell_types):
            info = signif_info.get(ev, {}).get(ct)
            if isinstance(info, dict):
                color = cmap(norm(float(np.clip(info['r'], -1, 1))))
                ec = 'black' if info.get('is_specific') else 'white'
                lw = 1.2    if info.get('is_specific') else 0.3
            else:
                color, ec, lw = '#eeeeee', 'white', 0.3

            ax.add_patch(Rectangle((j, i - 0.4), 1, 0.8,
                                    facecolor=color, edgecolor=ec,
                                    linewidth=lw, zorder=2))

            # dot for FDR < 0.01
            if isinstance(info, dict) and info.get('fdr', 1) < 0.01:
                ax.text(j + 0.5, i, '·', ha='center', va='center',
                        fontsize=9, color='white', zorder=3)

    ax.set_xlim(0, len(cell_types))
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xticks([j + 0.5 for j in range(len(cell_types))])
    ax.set_xticklabels(cell_types, rotation=40, ha='left', fontsize=5.5)
    ax.xaxis.set_tick_params(length=0)
    ax.tick_params(axis='x', pad=1)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=12)
    cbar.set_label('r', fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])


# ── main entry point ───────────────────────────────────────────────────────
def plot_cluster(cid, event_ids, exon_info, signif_info,
                 figsize=None, save_path=None, dpi=150):
    """
    Parameters
    ----------
    cid           : cluster id (for the title)
    event_ids     : members[cid]
    exon_info     : dict  event_id -> annotate_event record
    signif_info   : dict  event_id -> {'gene_name': str, cell_type: {'r','fdr','is_specific'}}
    figsize       : (width, height) in inches; None = auto
    save_path     : file path to save (e.g. 'cluster.pdf'); None = plt.show()
    dpi           : resolution for raster output
    """
    signif_set = set(signif_info.keys())
    rows       = _build_rows(event_ids, exon_info, signif_set)

    if not rows:
        print(f"cluster {cid}: no events found in exon_info")
        return None

    # cell types: all keys except 'gene_name', skip non-dict values defensively
    cell_types = sorted({
        ct
        for r in rows
        for ct, val in signif_info.get(r['event'], {}).items()
        if ct != 'gene_name' and isinstance(val, dict)
    })

    n      = len(rows)
    height = max(1.8, n * 0.38)

    if cell_types:
        n_ct         = len(cell_types)
        default_size = (7 + n_ct * 0.7, height)
        fig, (ax_l, ax_r) = plt.subplots(
            1, 2,
            figsize=figsize or default_size,
            gridspec_kw={'width_ratios': [4, max(1, n_ct)]}
        )
        _plot_celltype_panel(ax_r, rows, signif_info, cell_types)
        ax_r.set_title('cell-type associations\n(border = specific, · = FDR<0.01)',
                        fontsize=6, loc='left', pad=3)
    else:
        fig, ax_l = plt.subplots(figsize=figsize or (7, height))

    _plot_structure(ax_l, cid, rows, signif_set, signif_info)

    ax_l.legend(handles=[
        mpatches.Patch(facecolor='tomato',    label='significant (any ct)'),
        mpatches.Patch(facecolor='steelblue', label='not significant'),
    ], loc='upper right', fontsize=5.5, frameon=False)

    plt.tight_layout(pad=1.5)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"saved -> {save_path}")
    else:
        plt.show()

    return fig