import numpy as np
import matplotlib.pyplot as plt
from typing import Sequence, Protocol, List, TypeVar, Tuple

def _hamming_dist_matrix(seqs: list[str]) -> np.ndarray:
    """Normalized Hamming distance in [0,1]. Assumes all seqs same length."""
    if len(seqs) == 0:
        return np.zeros((0, 0), dtype=float)
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("All sequences must have the same length.")

    # Encode as uint8 array (N, L)
    A = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(len(seqs), L)
    D = np.empty((len(seqs), len(seqs)), dtype=float)
    for i in range(len(seqs)):
        # Vectorized mismatch count against row i
        mism = (A != A[i]).sum(axis=1)
        D[i] = mism / L
    return D

def _shannon_entropy_from_counts(counts: np.ndarray) -> float:
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

#TODO: This is a copy of the code in crossover.py. Should be moved to separate utils module and referenced.
def _parse_fasta(fasta_text: str) -> Tuple[str, str]:
    """
    Helper to parse a simple FASTA string into header and sequence.
    Assumes a single object with a single header."""
    if not fasta_text:
        return "", ""
    lines = fasta_text.strip().split('\n')
    if not lines:
        return "", ""
    header = lines[0] if lines[0].startswith('>') else ">"
    # The sequence could be on multiple lines or just one
    seq = "".join([line.strip() for line in lines if not line.startswith('>')])
    if not seq and not lines[0].startswith('>'):
        seq = lines[0].strip()
    return header, seq

def visualize_population_clusters(
    seqs: list[str],
    fitness: np.ndarray | None = None,
    gen: int | None = None,
    eps: float = 0.06,         # DBSCAN radius in normalized Hamming space
    min_samples: int = 3,      # DBSCAN core size
    embed: str = "pca",        # "pca" (default) or "umap" if you have umap-learn
    outpath: str | None = None,
    random_state: int = 0,
):
    """
    Produces a quick diversity snapshot:
      - 2D embedding (PCA or UMAP) colored by DBSCAN clusters
      - Species size bar chart
      - Summary metrics in title

    eps guideline (normalized Hamming):
      - 0.02 means ~2% of positions can differ
      - 0.06 means ~6% differ (often a good starting point for protein GA)
    """
    N = len(seqs)
    if N == 0:
        raise ValueError("Empty population.")

    fasta_seqs = seqs # Our inputs are actually FASTA-style strings
    seqs = []
    for f in fasta_seqs:
        _, seq = _parse_fasta(f)
        seqs.append(seq)

    D = _hamming_dist_matrix(seqs)  # (N,N)

    # --- clustering: DBSCAN on precomputed distances ---
    try:
        from sklearn.cluster import DBSCAN
    except ImportError as e:
        raise ImportError("Need scikit-learn for clustering (pip install scikit-learn).") from e

    labels = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed").fit_predict(D)
    # labels: -1 is "noise" (unassigned). Treat as its own color group(s) in the plot.

    # --- embedding ---
    if embed.lower() == "umap":
        try:
            import umap
        except ImportError:
            raise ImportError("embed='umap' requires umap-learn (pip install umap-learn).")
        # UMAP on distance matrix
        reducer = umap.UMAP(
            n_neighbors=min(15, max(5, N // 10)),
            min_dist=0.1,
            metric="precomputed",
            random_state=random_state,
        )
        XY = reducer.fit_transform(D)
    else:
        # PCA on classical MDS-style features via double-centering distances is more work.
        # For a quick stable view: do PCA on one-hot sequence features.
        # One-hot is big but fine at N<=500, L~few hundred.
        L = len(seqs[0])
        alphabet = sorted(set("".join(seqs)))
        a2i = {a: i for i, a in enumerate(alphabet)}
        K = len(alphabet)
        X = np.zeros((N, L * K), dtype=np.uint8)
        for n, s in enumerate(seqs):
            for j, aa in enumerate(s):
                X[n, j * K + a2i[aa]] = 1

        from sklearn.decomposition import PCA
        XY = PCA(n_components=2, random_state=random_state).fit_transform(X)

    # --- summary metrics ---
    iu = np.triu_indices(N, k=1)
    mean_pairwise = float(D[iu].mean()) if N > 1 else 0.0

    # cluster counts (exclude noise=-1 in the count list, but still report)
    unique, counts = np.unique(labels[labels >= 0], return_counts=True)
    n_clusters = int(len(unique))
    n_noise = int(np.sum(labels == -1))
    shannon = _shannon_entropy_from_counts(counts) if counts.size else 0.0
    max_cluster = int(counts.max()) if counts.size else 0

    # --- plotting ---
    fig = plt.figure(figsize=(10, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1.0], wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    # Color map by integer label; noise=-1 plotted in gray
    # We'll plot clusters first, noise last for visibility.
    if fitness is not None:
        fitness = np.asarray(fitness)
        # Filter nans:
        fitness[fitness == np.nan] = 0

    for lab in np.unique(labels[labels >= 0]):
        idx = np.where(labels == lab)[0]
        sizes = None
        if fitness is not None:
            # size scaled by fitness rank (safe + robust)
            r = fitness[idx].argsort().argsort()
            sizes = 30 + 70 * (r / max(1, (len(r) - 1)))
        ax.scatter(XY[idx, 0], XY[idx, 1], s=sizes if sizes is not None else 40, alpha=0.9, label=f"{lab}")
    
    if n_noise:
        idx = np.where(labels == -1)[0]
        ax.scatter(XY[idx, 0], XY[idx, 1], s=30, alpha=0.6, label="noise", color="gray")

    ax.set_xlabel("dim 1")
    ax.set_ylabel("dim 2")
    ax.set_title("Population embedding (colored by DBSCAN clusters)")

    # Species sizes bar chart (descending)
    if counts.size:
        sort = np.argsort(-counts)
        counts_sorted = counts[sort]
        axb.bar(np.arange(len(counts_sorted)), counts_sorted)
        axb.set_xticks([])
        axb.set_ylabel("count")
        axb.set_title("Species sizes")
    else:
        axb.text(0.5, 0.5, "No clusters\n(all noise)", ha="center", va="center")
        axb.set_axis_off()

    title = []
    if gen is not None:
        title.append(f"Gen {gen}")
    title.append(f"N={N}")
    title.append(f"clusters={n_clusters} (+noise={n_noise})")
    title.append(f"mean d={mean_pairwise:.3f}")
    title.append(f"H={shannon:.2f}")
    title.append(f"max cluster={max_cluster}")
    if fitness is not None:
        title.append(f"best fitness={float(np.max(fitness)):.3g}")
    fig.suptitle(" | ".join(title), y=1.02)

    ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()

    if outpath:
        fig.savefig(outpath, dpi=200, bbox_inches="tight")
        plt.close(fig)
    return fig, labels, {"mean_pairwise_dist": mean_pairwise, "shannon": shannon, "n_clusters": n_clusters, "n_noise": n_noise}