# ------------------------------------------------------------------
# 0.  Acknowledgment: This plotting code is implemented with the help of ChatGPT.
#     It was prompted to generate a hierarchical tree plot using matplotlib,
#     based on p-values from R output, with custom coloring and layout.
#     The code was then iteratively refined and adapted to our specific needs.
# ------------------------------------------------------------------

import argparse
import json
import re
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ------------------------------------------------------------------
# 1.  Utility helpers to build the tree dynamically from R output
# ------------------------------------------------------------------
level_box_width = {0: 14.0, 1: 12.0, 2: 10.0}
box_height = 1.6
leaf_marker_size = 900  # size for plt.scatter circles
cmap = cm.Blues_r
norm = mcolors.Normalize(vmin=0, vmax=1)


def parse_leaf_path(name: str):
    """
    Return the hierarchical path for a leaf based on its name.

    Expected pattern: Ages8to11_M_BMXBMI -> ["All Outcomes", "Ages 8-11", "M: 8-11", leaf]
    Falls back to targeted internal-node mappings when possible, otherwise
    ["All Outcomes", "Ungrouped", name]. Examples of targeted fallbacks:
      - Ages8to11          -> ["All Outcomes", "Ages 8-11"]
      - Ages8to11_m        -> ["All Outcomes", "Ages 8-11", "M: 8-11"]
      - Root / All Outcomes-> ["All Outcomes"]
    """
    match = re.match(r"^Ages(?P<rng>\d+to\d+)_([MF])_(?P<rest>.+)$", name)
    if match:
        rng_token = match.group("rng")
        rng_label = rng_token.replace("to", "-")
        sex = match.group(2)
        return ["All Outcomes", f"Ages {rng_label}", f"{sex}: {rng_label}", name], True

    no_suffix = re.match(r"^Ages(?P<rng>\d+to\d+)(?:_(?P<sex>[mMfF]))?$", name)
    if no_suffix:
        rng_token = no_suffix.group("rng")
        rng_label = rng_token.replace("to", "-")
        sex = no_suffix.group("sex")
        if sex:
            return ["All Outcomes", f"Ages {rng_label}", f"{sex.upper()}: {rng_label}"], False
        return ["All Outcomes", f"Ages {rng_label}"], False

    lowered = name.lower()
    if lowered in {"root", "all outcomes", "all_outcomes"}:
        return ["All Outcomes"], False

    return ["All Outcomes", "Ungrouped", name], True


def insert_path(tree: dict, path: list):
    """Insert a leaf path like ['All Outcomes', 'Ages 8-11', 'F: 8-11', 'leaf'] into the nested dict."""
    node = tree
    for i, segment in enumerate(path[:-1]):
        is_parent = i == len(path) - 2  # segment right before the leaf
        if is_parent:
            children = node.get(segment)
            if not isinstance(children, list):
                children = []
                node[segment] = children
            if path[-1] not in children:
                children.append(path[-1])
        else:
            if segment not in node or isinstance(node.get(segment), list):
                node[segment] = {}
            node = node[segment]
    return tree


def derive_tree_from_pvals(pvals: dict) -> dict:
    """Build a nested tree structure inferred from the named p-values."""
    tree = {}
    for name in pvals:
        path, is_leaf = parse_leaf_path(name)
        if is_leaf:  # only leaves determine the topology; internal nodes are implied
            insert_path(tree, path)
    sorted_tree = sort_tree({"All Outcomes": tree.get("All Outcomes", tree)})
    return sorted_tree


def normalize_pvals(pvals: dict) -> dict:
    """
    Map raw p-value keys coming from R to the node labels used in the plot.
    This aligns internal-node names like 'Ages8to11' with 'Ages 8-11'.
    """
    normalized = {}
    for name, val in pvals.items():
        path, is_leaf = parse_leaf_path(name)
        target = path[-1] if path else name
        normalized[target] = val
    return normalized


def _sort_key(label: str):
    """Custom sort key to keep ages increasing (8-11 before 12-17) and sexes grouped."""
    age_match = re.match(r"^Ages\s*(\d+)-(\d+)$", label)
    sex_match = re.match(r"^([MF]):\s*(\d+)-(\d+)$", label)
    if age_match:
        return (0, int(age_match.group(1)), int(age_match.group(2)), label)
    if sex_match:
        return (1, int(sex_match.group(2)), int(sex_match.group(3)), label)
    return (2, label)


def sort_tree(node):
    """Return a version of the tree with deterministic (sorted) ordering."""
    if isinstance(node, dict):
        ordered = {}
        for key in sorted(node, key=_sort_key):
            ordered[key] = sort_tree(node[key])
        return ordered
    if isinstance(node, list):
        return sorted(node, key=_sort_key)
    return node


def build_level_map(tree, level=0, mapping=None):
    """Return {node_name: depth} for every internal node in the tree."""
    mapping = mapping or {}
    if isinstance(tree, dict):
        for key, child in tree.items():
            mapping[key] = level
            build_level_map(child, level + 1, mapping)
    elif isinstance(tree, list):
        for child in tree:
            if isinstance(child, dict):
                build_level_map(child, level, mapping)
    return mapping


def positions_for_tree(tree, x0=0, y0=0, dx_internal=5, dx_leaf=3, dy=3):
    """Recursively assign x/y coordinates to every node/leaf."""
    if isinstance(tree, dict):
        pos = {}
        x = x0
        for k, sub in tree.items():
            subpos, x = positions_for_tree(sub, x, y0 - dy, dx_internal, dx_leaf, dy)
            cx_list = [pt[0] for pt in subpos.values()]
            pos[k] = (sum(cx_list) / len(cx_list), y0)
            pos.update(subpos)
        return pos, x + dx_internal
    elif isinstance(tree, list):
        pos = {leaf: (x0 + i * dx_leaf, y0) for i, leaf in enumerate(tree)}
        return pos, x0 + len(tree) * dx_leaf
    else:  # leaf string
        return {tree: (x0, y0)}, x0 + dx_internal


def pretty_leaf(label: str) -> str:
    """Drop the leading Ages*_Sex_ prefix for rotated labels."""
    bits = label.split("_")
    if len(bits) > 2 and bits[0].startswith("Ages"):
        return "_".join(bits[2:])
    return label


# ------------------------------------------------------------------
# 2.  Main plotting routine
# ------------------------------------------------------------------
def display_label(name: str, label_map: dict | None) -> str:
    """Return user-supplied label if available, otherwise the pretty leaf name."""
    if label_map and name in label_map:
        return label_map[name]
    return pretty_leaf(name)


def build_tree(gamma_label: str, tree: dict, pvals: dict, outfile: str, show_values: bool = False, label_map: dict | None = None):
    level_map = build_level_map(tree)
    positions, _ = positions_for_tree(tree)
    fig, ax = plt.subplots(figsize=(24, 14))

    def edges(node, parent=None):
        if isinstance(node, dict):
            for k, v in node.items():
                if parent:
                    px, py = positions[parent]
                    cx, cy = positions[k]
                    ax.plot(
                        [px, cx],
                        [py - box_height / 2, cy + box_height / 2],
                        color="gray",
                        lw=1.3,
                        zorder=0,
                    )
                edges(v, k)
        elif isinstance(node, list):
            for leaf in node:
                px, py = positions[parent]
                lx, ly = positions[leaf]
                ax.plot([px, lx], [py - box_height / 2, ly], color="gray", lw=1.2, zorder=0)

    edges(tree)

    for name, (x, y) in positions.items():
        pval = pvals.get(name, 1.0)
        if name in level_map:
            width = level_box_width.get(level_map[name], 10.0)
            face = cmap(norm(pval))
            r, g, b, _ = mcolors.to_rgba(face)
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            box = FancyBboxPatch(
                (x - width / 2, y - box_height / 2),
                width,
                box_height,
                boxstyle="round,pad=0.3",
                fc=face,
                ec="black",
                lw=1.6,
                zorder=2,
            )
            ax.add_patch(box)
            ax.text(
                x,
                y,
                name,
                ha="center",
                va="center",
                fontsize=22,
                weight="bold",
                color="white" if luminance < 0.5 else "black",
                zorder=3,
            )
            if show_values:
                ax.text(
                    x,
                    y - box_height / 2 - 0.35,
                    f"p={pval:.4g}",
                    ha="center",
                    va="top",
                    fontsize=12,
                    color="black",
                    zorder=3,
                )
        else:
            ax.scatter(x, y, s=leaf_marker_size, c=[cmap(norm(pval))], edgecolors="black", zorder=2)
            if pval < 0.05:
                ax.text(
                    x,
                    y - 1.0,
                    display_label(name, label_map),
                    ha="center",
                    va="top",
                    rotation=90,
                    fontsize=18,
                    weight="bold",
                    color="#000000",
                    zorder=3,
                )
            if show_values:
                ax.text(
                    x,
                    y + 0.3,
                    f"{pval:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    color="#000000",
                    rotation=90,
                    zorder=3,
                )

    ax.set_title(fr"$\Gamma = {gamma_label}$", fontsize=55, y=0.875, weight="bold")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.046, pad=0.08)
    cbar.set_label("p-value", fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    xs, ys = zip(*positions.values())
    ax.set_xlim(min(xs) - 7, max(xs) + 7)
    ax.set_ylim(min(ys) - 5, max(ys) + 5)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {outfile}")


# ------------------------------------------------------------------
# 3.  Loading R output and CLI entry point
# ------------------------------------------------------------------
def load_r_output(json_path: Path) -> dict:
    """Load list-like JSON from R and return {gamma_value: pvals_dict}."""
    with json_path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    results = {}
    for entry in payload:
        gamma_val = entry["gamma"]
        results[gamma_val] = entry["pvalues"]
    return results


def load_labels(path: Path | None) -> dict:
    """Load optional JSON mapping from raw variable names to display labels."""
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def main():
    parser = argparse.ArgumentParser(description="Render outcome tree from R-generated p-values.")
    parser.add_argument(
        "--gamma",
        type=float,
        help="Gamma value to plot (use --gammas for multiple).",
    )
    parser.add_argument(
        "--gammas",
        type=float,
        nargs="+",
        help="List of Gamma values to plot (overrides --gamma if provided).",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Path to JSON produced by R (list of {gamma, pvalues}). If omitted, will try r_output/pvalues_gamma<gamma>.json then r_output/pvalues.json.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Where to write the plot. If omitted, saves to <prefix>outcomes_tree_gamma<gamma>.pdf.",
    )
    parser.add_argument(
        "--show-values",
        action="store_true",
        help="Overlay p-values on the plot for audit/debugging.",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Optional JSON mapping of raw variable names to display labels.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="",
        help='Optional prefix applied to default input/output filenames (e.g., "unadj_" or "bonf_").',
    )
    args = parser.parse_args()

    gamma_list = args.gammas or ([args.gamma] if args.gamma is not None else None)
    if not gamma_list:
        raise SystemExit("Provide --gamma <val> or --gammas <v1> <v2> ...")

    try:
        gamma_values = [float(g) for g in gamma_list]
    except (TypeError, ValueError):
        raise SystemExit(f"All gamma values must be numeric; received {gamma_list!r}")

    label_map = load_labels(args.labels or Path("r_output/labels.json"))

    for gamma_val in gamma_values:
        gamma_label = f"{gamma_val:.2f}"
        input_path = args.input
        if input_path is None:
            guess = Path(f"r_output/{args.prefix}pvalues_gamma{gamma_label}.json")
            input_path = guess if guess.exists() else Path(f"r_output/{args.prefix}pvalues.json")

        if args.output:
            output_path = args.output
            if len(gamma_values) > 1:
                output_path = output_path.with_name(
                    f"{output_path.stem}_{gamma_label}{output_path.suffix}"
                )
        else:
            output_path = Path(f"{args.prefix}outcomes_tree_gamma{gamma_label}.pdf")

        print(f"Reading p-values from: {input_path} for Gamma={gamma_label}")
        results = load_r_output(input_path)
        gamma_key = None
        if gamma_val in results:
            gamma_key = gamma_val
        else:
            # tolerate tiny float differences or stringified keys
            for k in results:
                if isinstance(k, (int, float)) and isinstance(gamma_val, (int, float)):
                    try:
                        if abs(k - gamma_val) < 1e-9:
                            gamma_key = k
                            break
                    except TypeError:
                        pass
                if str(k) == str(gamma_val):
                    gamma_key = k
                    break
        if gamma_key is None:
            available = ", ".join(str(g) for g in sorted(results.keys()))
            print(f"Gamma {gamma_val} not found in input {input_path}. Available: {available}")
            continue

        pvals = results[gamma_key]
        tree = derive_tree_from_pvals(pvals)
        aligned_pvals = normalize_pvals(pvals)
        build_tree(
            f"{float(gamma_key):.2f}" if isinstance(gamma_key, (int, float)) else str(gamma_key),
            tree,
            aligned_pvals,
            str(output_path),
            show_values=args.show_values,
            label_map=label_map,
        )
        print(f"Saved plot to: {output_path}")


if __name__ == "__main__":
    main()
