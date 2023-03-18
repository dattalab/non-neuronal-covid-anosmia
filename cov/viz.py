import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


def plot_heatmap(
    adata, genes_of_interest, cov_genes, order, clust="cluster_name"
):
    """
    Plot a heatmap of normalized expression for a set of genes.

    Args:
        adata (anndata.AnnData)
        genes_of_interest (list): list of cell-type markers
        cov_genes (list): list of sars-related genes to plot
        subset_order (list): ordered subset of cell types to use
        clust (str, optional): column in adata.obs with cluster labels. Defaults to "cluster_name".

    Returns:
        df of shape (clusters x genes) with expression normalized across clusters for each gene.
    """
    n_goi = len(genes_of_interest)
    n_cov = len(cov_genes)
    either_gene = genes_of_interest + cov_genes
    df_genes = sc.get.obs_df(adata, either_gene)
    df_mean = df_genes.join(adata.obs[clust]).groupby(clust).mean()
    mean_df_sub = df_mean.loc[order]
    mean_df_sub_norm = mean_df_sub / mean_df_sub.max(axis=0)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.3, 4.8),
        gridspec_kw={"width_ratios": [n_cov, n_goi]},
        sharey=True,
    )
    cb_ax = fig.add_axes([0.91, 0.12, 0.015, 0.75])
    ax1, ax2 = axes

    sns.heatmap(mean_df_sub_norm[cov_genes], cmap="Reds", ax=ax1, cbar=False)
    hm = sns.heatmap(
        mean_df_sub_norm[genes_of_interest],
        cmap="Reds",
        ax=ax2,
        cbar_ax=cb_ax,
        cbar_kws={"label": ("Normalized expression")},
    )
    fig.subplots_adjust(wspace=0.02)
    ax2.tick_params(left=False)
    ax2.set_ylabel(None)
    ax1.set_ylabel("Cluster")
    return mean_df_sub_norm


def make_bar_plot(
    adata,
    clusts,
    labs,
    viral_genes=["ACE2", "TMPRSS2"],
    clust="cluster_name",
    ylims=(4, 60),
    mouse=False,
):
    """
    Make plot of the percentage of cell types expressing viral genes.

    Args:
        adata (anndata.AnnData)
        clusts (tuple): tuple of (olf, resp) cluster names to use
        labs (tuple): tuple (olf, resp) xticklabels for each cluster
        viral_genes (list, optional): viral genes to plot. Defaults to ["ACE2", "TMPRSS2"].
        clust (str, optional): column in adata.obs with cluster labels. Defaults to "cluster_name".
        ylims (tuple, optional): ylims for each viral gene
        mouse (bool, optional): whether it's mouse data (converts gene names)

    Returns:
        df: dataframe of the percent of cells expressing each viral gene for each cluster
    """
    if mouse:
        viral_genes = [v.title() for v in viral_genes]
    AC, TM = viral_genes

    pct_exp = (sc.get.obs_df(adata, viral_genes) > 0).join(adata.obs[clust]).groupby(
        clust
    ).mean() * 100
    df_to_plot = pct_exp
    cls = plt.cm.Set1.colors

    non_resp, resp = clusts
    n_non = len(non_resp)
    n_resp = len(resp)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8, 4),
        gridspec_kw={"width_ratios": [n_non, n_resp]},
        sharey=False,
    )

    for ax, subset, title, labels in zip(
        axes.flatten(), (non_resp, resp), ["Olfactory", "Respiratory"], labs
    ):
        axb = ax.twinx()
        width = 0.4
        df_to_plot[AC].loc[subset].plot(
            kind="bar", color=cls[0], ax=ax, width=width, position=1
        )
        df_to_plot.loc[subset][TM].plot(
            kind="bar", color=cls[1], ax=axb, width=width, position=0
        )
        ax.set_xticklabels(labels, rotation=90)
        ax.text(
            0,
            1.03,
            f"% {AC}+",
            horizontalalignment="center",
            verticalalignment="bottom",
            transform=ax.transAxes,
            fontdict={"color": cls[0]},
        )
        axb.text(
            1,
            1.03,
            f"% {TM}+",
            horizontalalignment="center",
            verticalalignment="bottom",
            transform=ax.transAxes,
            fontdict={"color": cls[1]},
        )
        lims = list(ax.get_xlim())
        lims[0] = -0.5
        ax.set_xlim(lims)
        axb.set_ylim(0, ylims[1])
        ax.set_ylim(0, ylims[0])
        ax.set_title(title, y=1.08)
        ax.tick_params(axis="y", colors=cls[0])
        axb.tick_params(axis="y", colors=cls[1])
        axb.spines["left"].set_color(cls[0])
        axb.spines["right"].set_color(cls[1])
        ax.set_xlabel(None)
    sns.despine(right=False)
    plt.subplots_adjust(wspace=0.5)
    return df_to_plot


def make_color_dict(gen_clust_order):
    cell_colors = {
        "Olf. HBC": (0.0902, 0.7451, 0.81176),
        "OSN": (0.49804, 0.49804, 0.49804),
        "SUS": (0.8902, 0.46667, 0.76078),
        "Bowman's gland": (0.76863, 0.61176, 0.58039),
        "MV Brush-like": (0.54902, 0.33725, 0.29412),
        "MV Ionocyte-like": (0.77255, 0.6902, 0.83529),
        "OEC": (0.58039, 0.40392, 0.74118),
        "Resp. HBC": (1.0, 0.59608, 0.58824),
        "Resp. secretory": (0.59608, 0.87451, 0.54118),
        "Resp. ciliated": (0.17255, 0.62745, 0.17255),
        "Erythrocyte": (0.68235, 0.78039, 0.9098),
    }
    iter_cmap = iter(list(plt.cm.Dark2.colors) + list(plt.cm.Set3.colors[2:]))
    col_dict = {}
    for c in gen_clust_order:
        if c in cell_colors:
            col_dict[c] = cell_colors[c]
        else:
            col_dict[c] = next(iter_cmap)

    fix_col = lambda t: [x / 256 for x in t]
    col_dict["T cell"] = fix_col((253, 192, 134))
    col_dict["Cycling Resp. HBC"] = fix_col((240, 228, 66))
    col_dict["MV Brush-like"] = fix_col((227, 26, 28))
    col_dict["Circulating myeloid"] = fix_col((44, 69, 125))
    col_dict["Resp. HBC"] = fix_col((163, 21, 38))
    col_dict["B cell"] = fix_col((255, 179, 208))
    return col_dict


def hbc_colors():
    color_dict = {
        "HBC": (0.0902, 0.7451, 0.81176),
        "HBC*": (1.0, 0.49804, 0.0549),
        "GBC": (0.85882, 0.85882, 0.55294),
        "INP": (0.73725, 0.74118, 0.13333),
        "iOSN": (0.78039, 0.78039, 0.78039),
        "mOSN": (0.49804, 0.49804, 0.49804),
        "SUS": (0.8902, 0.46667, 0.76078),
        "MV Ionocyte-like": (0.58039, 0.40392, 0.74118),
        "Resp. HBC": (1.0, 0.59608, 0.58824),
        "Resp. Secretory": (0.17255, 0.62745, 0.17255),
    }
    return color_dict


def add_cbar(fig, scatter, pos=None):
    """
    Add a vertical colorbar with updated ticks

    Args:
        fig: mpl figure
        scatter: dataset to add colorbar to (output of ax.scatter)
        pos (list, optional): coordinates for colorbar. Defaults to None.

    Returns:
        cbar
    """
    cb_ax = fig.add_axes(pos)
    cbar = fig.colorbar(scatter, cax=cb_ax)
    cbar.set_label("Normalized expression", labelpad=-10)
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")
    cbar.set_ticks(cbar.mappable.get_clim())
    cbar.set_ticklabels([0, "max"])
    return cbar