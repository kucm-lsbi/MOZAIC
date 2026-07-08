import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------- Style ---------------------------------------- #

PLOT_STYLE = {
    "font.family": "DejaVu Sans",
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",    
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "font.weight": "bold",            
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.alpha": 0.30,
    "axes.axisbelow": True,
    "figure.dpi": 110,
}

METRIC_COLORS = {
    "qed":  "#d62728",   
    "sa":   "#2ca02c",   
    "vina": "#1f77b4",   
}

def categorical_palette(n):
    base = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#393b79",
    ]
    return [base[i % len(base)] for i in range(n)]


# ----------------------------- Pareto helpers ------------------------------- #

def pareto_mask(scores_to_max: np.ndarray) -> np.ndarray:
    """
    Return a boolean mask of non-dominated points.
    scores_to_max: (n_points, n_obj) — all objectives oriented as 'higher is better'.
    """
    n = scores_to_max.shape[0]
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        if not mask[i]:
            continue
        diff = scores_to_max - scores_to_max[i]
        dominated_by = np.all(diff >= 0, axis=1) & np.any(diff > 0, axis=1)
        if np.any(dominated_by):
            mask[i] = False
    return mask


def hypervolume_2d(points: np.ndarray, ref: np.ndarray) -> float:

    if points.size == 0:
        return 0.0
    pts = points[points[:, 0].argsort()[::-1]]
    hv = 0.0
    prev_y = ref[1]
    for x, y in pts:
        if y > prev_y:
            hv += (x - ref[0]) * (y - prev_y)
            prev_y = y
    return float(hv)


# ----------------------------- Tracker -------------------------------------- #

class CSATracker:
    """
    Same `track_bank(bank, round_num, davg=None)` signature.
    Produces 4 publication-clean figures with simplified titles and breathing layouts.
    """

    def __init__(self, job_dir: str):
        self.job_dir = job_dir
        os.makedirs(job_dir, exist_ok=True)
        self.history = []
        self.snapshots = []
        self.round_boundaries = []
        self.step_counter = -1

    # ---------- ingestion ---------- #

    @staticmethod
    def _clean_rxn_name(pos):
        name = re.sub(r"^\d+_", "", str(pos))
        name = re.sub(r"_\d+$", "", name).replace("_", " ").strip()
        return name or "unknown"

    @staticmethod
    def _scale_to_0_1(x):
        x = np.asarray(x, dtype=float)
        lo, hi = np.nanmin(x), np.nanmax(x)
        if not np.isfinite(lo) or not np.isfinite(hi) or abs(hi - lo) < 1e-12:
            return np.zeros_like(x, dtype=float)
        return (x - lo) / (hi - lo)

    def track_bank(self, bank, round_num, davg=None):

        self.step_counter += 1
        if round_num > 0 and (not self.round_boundaries
                              or self.round_boundaries[-1]["round"] != round_num):
            self.round_boundaries.append({"round": round_num, "step": self.step_counter})

        vinas = np.array([m.get("scores", {}).get("affinity", np.nan) for m in bank], dtype=float)
        qeds  = np.array([m.get("scores", {}).get("qed",      np.nan) for m in bank], dtype=float)
        sa_raw = np.array([m.get("scores", {}).get("sa", np.nan) for m in bank], dtype=float)
        sas = -sa_raw

        self.snapshots.append({"step": self.step_counter, "round": round_num,
                               "qed": qeds, "sa": sas, "sa_raw": sa_raw,
                               "vina": vinas})

        top_k = max(1, int(len(bank) * 0.10))
        rank = (-vinas) * 0.50 + qeds * 2.00 + sas * 2.00
        rank[np.isnan(rank)] = -np.inf
        top_idx = set(np.argsort(-rank)[:top_k])

        rxn_counts, rxn_top = defaultdict(int), defaultdict(int)
        rxn_scores = defaultdict(lambda: {"qed": [], "sa": [], "vina": []})

        for i, m in enumerate(bank):
            mol_qed, mol_sa, mol_vina = qeds[i], sas[i], vinas[i]
            for r in m.get("rxn_history", []):
                if r.get("selected_fragment") == "skip":
                    continue
                name = self._clean_rxn_name(r.get("rxn_position", "unknown"))
                rxn_counts[name] += 1
                if i in top_idx:
                    rxn_top[name] += 1
                if np.isfinite(mol_qed):  rxn_scores[name]["qed"].append(mol_qed)
                if np.isfinite(mol_sa):   rxn_scores[name]["sa"].append(mol_sa)
                if np.isfinite(mol_vina): rxn_scores[name]["vina"].append(mol_vina)

        sa_abs_min  = -np.nanmax(sas)
        sa_abs_max  = -np.nanmin(sas)
        sa_abs_mean = -np.nanmean(sas)

        self.history.append({
            "step": self.step_counter, "round": round_num,
            "qed_min": np.nanmin(qeds), "qed_mean": np.nanmean(qeds),
            "qed_max": np.nanmax(qeds), "qed_std": np.nanstd(qeds),
            "sa_min": np.nanmin(sas),   "sa_mean":  np.nanmean(sas),
            "sa_max": np.nanmax(sas),   "sa_std":   np.nanstd(sas),
            "sa_abs_min":  sa_abs_min,  "sa_abs_mean": sa_abs_mean,
            "sa_abs_max":  sa_abs_max,
            "vina_min": np.nanmin(vinas), "vina_mean": np.nanmean(vinas),
            "vina_max": np.nanmax(vinas), "vina_std":  np.nanstd(vinas),
            "davg": float(davg) if davg is not None else np.nan,
            "rxn_counts": dict(rxn_counts), "rxn_top_counts": dict(rxn_top),
            "rxn_scores": {k: {kk: list(vv) for kk, vv in v.items()}
                           for k, v in rxn_scores.items()},
        })

    # ---------- frames ---------- #

    def _summary_frames(self):
        df = pd.DataFrame([{k: v for k, v in r.items()
                            if k not in ("rxn_counts", "rxn_top_counts", "rxn_scores")}
                           for r in self.history])
        steps = df["step"].values

        rxn_data = pd.DataFrame([r["rxn_counts"] for r in self.history]).fillna(0)
        rxn_data.index = steps
        rxn_data = rxn_data.loc[:, rxn_data.sum(axis=0).sort_values(ascending=False).index]
        rxn_shares = rxn_data.div(rxn_data.sum(axis=1).replace(0, 1), axis=0)

        top_data = pd.DataFrame([r["rxn_top_counts"] for r in self.history]).fillna(0)
        top_data.index = steps
        top_data = top_data.reindex(columns=rxn_data.columns, fill_value=0)

        all_total, top_total = rxn_data.sum(axis=0), top_data.sum(axis=0)
        all_share = all_total / max(all_total.sum(), 1)
        top_share = top_total / max(top_total.sum(), 1)
        enrichment = (top_share / all_share.replace(0, np.nan)).replace(
            [np.inf, -np.inf], np.nan).fillna(0)

        contribution = pd.DataFrame({
            "all_count": all_total, "top_count": top_total,
            "all_share": all_share, "top_share": top_share,
            "top_enrichment": enrichment,
        }).fillna(0)

        rxn_delta = pd.DataFrame({
            "initial_share": rxn_shares.iloc[0],
            "final_share":   rxn_shares.iloc[-1],
            "delta":         rxn_shares.iloc[-1] - rxn_shares.iloc[0],
        }).sort_values("delta", ascending=True)

        return df, rxn_shares, contribution, rxn_delta

    def _reaction_coupling_frame(self):
        pool = {"qed": [], "sa": [], "vina": []}
        for h in self.history:
            for rs in h["rxn_scores"].values():
                pool["qed"].extend(rs["qed"])
                pool["sa"].extend(rs["sa"])
                pool["vina"].extend(rs["vina"])
        ref = {k: (np.mean(v) if v else 0.0,
                   np.std(v)  if v and np.std(v) > 1e-9 else 1.0)
               for k, v in pool.items()}

        agg = defaultdict(lambda: {"qed": [], "sa": [], "vina": []})
        for h in self.history:
            for r, rs in h["rxn_scores"].items():
                agg[r]["qed"].extend(rs["qed"])
                agg[r]["sa"].extend(rs["sa"])
                agg[r]["vina"].extend(rs["vina"])

        rows = []
        for r, rs in agg.items():
            n = max(len(rs["qed"]), len(rs["sa"]), len(rs["vina"]))
            if n == 0:
                continue
            mq = np.mean(rs["qed"]) if rs["qed"] else np.nan
            ms = np.mean(rs["sa"])  if rs["sa"]  else np.nan
            mv = np.mean(rs["vina"]) if rs["vina"] else np.nan
            rows.append({
                "reaction": r, "n": n,
                "mean_qed": mq, "mean_sa": ms, "mean_vina": mv,
                "qed_z":  (mq - ref["qed"][0])  / ref["qed"][1],
                "sa_z":   (ms - ref["sa"][0])   / ref["sa"][1],
                "vina_z": -(mv - ref["vina"][0]) / ref["vina"][1],
            })
        out = pd.DataFrame(rows).set_index("reaction")
        out["overall_z"] = out[["qed_z", "sa_z", "vina_z"]].sum(axis=1)
        out = out.sort_values("overall_z", ascending=True)
        return out

    def _objective_frames(self, df):
        steps = df["step"].values
        raw_delta = pd.DataFrame({
            "step": steps,
            "QED": df["qed_mean"].diff().fillna(0),
            "SA":  df["sa_mean"].diff().fillna(0),
            "Vina": -df["vina_mean"].diff().fillna(0),
        })
        scaled = pd.DataFrame({
            "step": steps,
            "QED":  self._scale_to_0_1(df["qed_mean"]),
            "SA":   self._scale_to_0_1(df["sa_mean"]),
            "Vina": self._scale_to_0_1(-df["vina_mean"]),
        })
        scaled_delta = scaled[["QED", "SA", "Vina"]].diff().fillna(0)

        drivers, colors = [], []
        for _, row in scaled_delta.iterrows():
            v = row[["QED", "SA", "Vina"]]
            if v.max() <= 0:
                drivers.append("None"); colors.append("#cccccc")
            else:
                d = v.idxmax(); drivers.append(d)
                colors.append(METRIC_COLORS["qed"] if d == "QED"
                              else METRIC_COLORS["sa"] if d == "SA"
                              else METRIC_COLORS["vina"])
        return raw_delta, scaled, pd.DataFrame({"step": steps, "driver": drivers, "color": colors})

    def _pareto_frames(self):
        all_qed  = np.concatenate([s["qed"]  for s in self.snapshots])
        all_sa   = np.concatenate([s["sa"]   for s in self.snapshots])
        all_vina = np.concatenate([-s["vina"] for s in self.snapshots])

        ref_qed  = np.nanmin(all_qed)  - 0.02
        ref_sa   = np.nanmin(all_sa)   - 0.02
        ref_vina = np.nanmin(all_vina) - 0.10

        qed_lo,  qed_hi  = np.nanmin(all_qed),  np.nanmax(all_qed)
        sa_lo,   sa_hi   = np.nanmin(all_sa),   np.nanmax(all_sa)
        vina_lo, vina_hi = np.nanmin(all_vina), np.nanmax(all_vina)
        def _safe_range(lo, hi):
            return (hi - lo) if (hi - lo) > 1e-9 else 1.0
        qed_rng  = _safe_range(qed_lo, qed_hi)
        sa_rng   = _safe_range(sa_lo,  sa_hi)
        vina_rng = _safe_range(vina_lo, vina_hi)

        ref_norm = -0.02

        rows = []
        for snap in self.snapshots:
            qed, sa, vina_max = snap["qed"], snap["sa"], -snap["vina"]
            scores3d = np.column_stack([qed, sa, vina_max])
            mask3d = pareto_mask(scores3d)

            qed_n  = (qed       - qed_lo)  / qed_rng
            sa_n   = (sa        - sa_lo)   / sa_rng
            vina_n = (vina_max  - vina_lo) / vina_rng

            pairs_raw = {
                "QED-Vina": (qed, vina_max, ref_qed, ref_vina),
                "SA-Vina":  (sa,  vina_max, ref_sa,  ref_vina),
                "QED-SA":   (qed, sa,       ref_qed, ref_sa),
            }
            pairs_norm = {
                "QED-Vina": (qed_n, vina_n),
                "SA-Vina":  (sa_n,  vina_n),
                "QED-SA":   (qed_n, sa_n),
            }

            hvs_raw, hvs_norm, fronts = {}, {}, {}
            for k in pairs_raw:
                a, b, ra, rb = pairs_raw[k]
                pts = np.column_stack([a, b])
                m = pareto_mask(pts)
                fronts[k] = pts[m]
                hvs_raw[k] = hypervolume_2d(pts[m], np.array([ra, rb]))

                an, bn = pairs_norm[k]
                pts_n = np.column_stack([an, bn])
                m_n = pareto_mask(pts_n)
                hvs_norm[k] = hypervolume_2d(pts_n[m_n], np.array([ref_norm, ref_norm]))

            rows.append({
                "step": snap["step"], "round": snap["round"],
                "n_pareto_3d": int(mask3d.sum()),
                "hv_QED_Vina_raw": hvs_raw["QED-Vina"],
                "hv_SA_Vina_raw":  hvs_raw["SA-Vina"],
                "hv_QED_SA_raw":   hvs_raw["QED-SA"],
                "hv_QED_Vina": hvs_norm["QED-Vina"],
                "hv_SA_Vina":  hvs_norm["SA-Vina"],
                "hv_QED_SA":   hvs_norm["QED-SA"],
                "fronts": fronts,
            })
        return rows

    # ---------- plotting ---------- #

    def _round_lines(self, ax):
        for rb in self.round_boundaries:
            ax.axvline(rb["step"], color="0.25", ls=":", lw=1.4, alpha=0.85)

    def _round_labels(self, ax, y=None):
        for rb in self.round_boundaries:
            ax.text(rb["step"], 1.02, f"R{rb['round']}",
                    transform=ax.get_xaxis_transform(),
                    rotation=0, va="bottom", ha="center",
                    fontweight="bold", fontsize=10, color="0.15",
                    clip_on=False)

    # ============ Objective Dynamics ============ #
    def fig_objective_dynamics(self, df, raw_delta, scaled, driver_df, save_path):
        steps = df["step"].values
        red, grn, blu = METRIC_COLORS["qed"], METRIC_COLORS["sa"], METRIC_COLORS["vina"]

        has_diversity = "davg" in df.columns and df["davg"].notna().any()

        if has_diversity:
            fig = plt.figure(figsize=(11, 13), constrained_layout=False)
            fig.subplots_adjust(left=0.10, right=0.97, top=0.96, bottom=0.05,
                                hspace=0.85)
            gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.0, 0.18, 0.7])
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = fig.add_subplot(gs[3])
        else:
            fig = plt.figure(figsize=(11, 11), constrained_layout=False)
            fig.subplots_adjust(left=0.10, right=0.97, top=0.96, bottom=0.06,
                                hspace=0.85)
            gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.18])
            ax1 = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
            ax3 = fig.add_subplot(gs[2])
            ax4 = None

        x_lo, x_hi = steps[0] - 1.8, steps[-1] + 0.5

        qed_raw = df["qed_mean"].values
        sa_internal = df["sa_mean"].values
        vina_raw = df["vina_mean"].values

        qed_lo, qed_hi   = float(np.nanmin(qed_raw)),  float(np.nanmax(qed_raw))
        sa_worst = -float(np.nanmin(sa_internal))
        sa_best  = -float(np.nanmax(sa_internal))
        vina_worst = float(np.nanmax(vina_raw))
        vina_best  = float(np.nanmin(vina_raw))

        qed_label  = f"QED\n[{qed_lo:.2f} → {qed_hi:.2f}]"
        sa_label   = f"SA\n[{sa_worst:.2f} → {sa_best:.2f}]"
        vina_label = f"Vina\n[{vina_worst:.2f} → {vina_best:.2f}]"

        ax1.plot(steps, scaled["QED"],  color=red, marker="o", lw=2.2, label=qed_label)
        ax1.plot(steps, scaled["SA"],   color=grn, marker="s", lw=2.2, label=sa_label)
        ax1.plot(steps, scaled["Vina"], color=blu, marker="^", lw=2.2, label=vina_label)
        ax1.set_title("(a) Run-wide rescaled bank mean progress", pad=22)
        ax1.set_ylabel("Normalized to run range")
        ax1.set_ylim(-0.05, 1.08)
        ax1.set_xticks(steps); ax1.set_xlim(x_lo, x_hi)
        ax1.legend(ncol=1, loc="upper left", frameon=False, labelspacing=1.6)
        self._round_lines(ax1)

        # (b) Per-update rescaled Δ
        scaled_delta = scaled[["QED", "SA", "Vina"]].diff().fillna(0)
        w = 0.27
        ax2.bar(steps - w, scaled_delta["QED"],  width=w, color=red, alpha=0.85, label="ΔQED")
        ax2.bar(steps,     scaled_delta["SA"],   width=w, color=grn, alpha=0.85, label="ΔSA")
        ax2.bar(steps + w, scaled_delta["Vina"], width=w, color=blu, alpha=0.85, label="ΔVina")
        ax2.axhline(0, color="black", lw=0.8)
        ax2.set_title("(b) Per-update rescaled improvement", pad=12)
        ax2.set_ylabel("Δ (rescaled)")
        ax2.set_xticks(steps); ax2.set_xlim(x_lo, x_hi)
        ax2.legend(ncol=1, loc="upper left", frameon=False)
        self._round_lines(ax2)

        # (c) Driver track
        ax3.bar(steps, np.ones(len(steps)), color=driver_df["color"].values,
                width=0.85, edgecolor="white", lw=0.5)
        full_label = {"QED": "QED", "SA": "SA", "Vina": "Vina", "None": "—"}
        for s, d in zip(driver_df["step"], driver_df["driver"]):
            ax3.text(s, 0.5, full_label[d], ha="center", va="center",
                     color="white" if d != "None" else "0.4",
                     fontweight="bold", fontsize=9)
        ax3.set_yticks([])
        ax3.set_xticks(steps); ax3.set_xlim(x_lo, x_hi)
        ax3.set_title("(c) Dominant improving objective per update", pad=12)
        ax3.set_ylim(0, 1)
        ax3.grid(False)
        self._round_lines(ax3)

        legend_handles = [
            Patch(facecolor=red, label="QED"),
            Patch(facecolor=grn, label="SA"),
            Patch(facecolor=blu, label="Vina"),
            Patch(facecolor="#cccccc", label="No gain"),
        ]
        ax3.legend(handles=legend_handles, ncol=4, loc="upper center",
                   bbox_to_anchor=(0.5, -0.7), frameon=False)

        if not has_diversity:
            ax3.set_xlabel("Update step")

        # (d) Bank distance average
        if has_diversity:
            davg_vals = df["davg"].values
            initial_davg = davg_vals[0]
            final_davg = davg_vals[-1]

            ax4.plot(steps, davg_vals, color="#6a3d9a", marker="D",
                     lw=2.2, markersize=6, label="Davg", zorder=2)

            init_color = "#1f9e3d"
            fin_color  = "#d62728"
            ax4.scatter([steps[0]], [initial_davg], s=130, color=init_color,
                        edgecolor="white", lw=1.2, zorder=4,
                        label=f"initial Davg\n[{initial_davg:.3f}]")
            ax4.scatter([steps[-1]], [final_davg], s=130, color=fin_color,
                        edgecolor="white", lw=1.2, zorder=4,
                        label=f"final Davg\n[{final_davg:.3f}]")

            ax4.set_xlabel("Update step")
            ax4.set_ylabel("Davg")
            ax4.set_title("(d) Bank distance average", pad=12)
            ax4.set_xticks(steps); ax4.set_xlim(x_lo, x_hi)
            dmin, dmax = float(np.nanmin(davg_vals)), float(np.nanmax(davg_vals))
            pad = max(0.005, (dmax - dmin) * 0.15)
            ax4.set_ylim(dmin - pad, dmax + pad)
            ax4.legend(loc="center left", frameon=False, fontsize=9,
                       labelspacing=1.6)
            self._round_lines(ax4)

        self._round_labels(ax1)

        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ============ Pareto Dynamics ============ #
    def fig_pareto(self, pareto_rows, save_path):
        steps = [r["step"] for r in pareto_rows]
        n_steps = len(steps)
        idx_init, idx_fin = 0, n_steps - 1

        fig = plt.figure(figsize=(13, 11), constrained_layout=False)
        fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.06,
                            wspace=0.32, hspace=0.55)
        gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 0.85])

        pair_keys = ["QED-Vina", "SA-Vina", "QED-SA"]
        pair_labels = {
            "QED-Vina": ("QED", "−Vina (kcal/mol)"),
            "SA-Vina":  ("SA (lower = easier)", "−Vina (kcal/mol)"),
            "QED-SA":   ("QED", "SA (lower = easier)"),
        }

        def get_xy_for_pareto(snap_idx, key):
            snap = self.snapshots[snap_idx]
            if key == "QED-Vina": return snap["qed"], -snap["vina"]
            if key == "SA-Vina":  return snap["sa"],  -snap["vina"]
            return snap["qed"], snap["sa"]

        def get_xy_for_display(snap_idx, key):
            snap = self.snapshots[snap_idx]
            sa_abs = snap.get("sa_raw", -snap["sa"])
            if key == "QED-Vina": return snap["qed"], -snap["vina"]
            if key == "SA-Vina":  return sa_abs,      -snap["vina"]
            return snap["qed"], sa_abs

        shared_lims = {}
        for key in pair_keys:
            xi, yi = get_xy_for_display(idx_init, key)
            xf, yf = get_xy_for_display(idx_fin,  key)
            xs = np.concatenate([xi, xf]); ys = np.concatenate([yi, yf])
            pad_x = (xs.max() - xs.min()) * 0.05
            pad_y = (ys.max() - ys.min()) * 0.05
            shared_lims[key] = ((xs.min() - pad_x, xs.max() + pad_x),
                                (ys.min() - pad_y, ys.max() + pad_y))

        # Rows 0–1: 2D pairwise scatter
        for col, key in enumerate(pair_keys):
            xlim, ylim = shared_lims[key]
            for row, snap_idx, tag in [(0, idx_init, "Initial"),
                                        (1, idx_fin,  "Final")]:
                ax = fig.add_subplot(gs[row, col])
                a_pf, b_pf = get_xy_for_pareto(snap_idx, key)
                a_disp, b_disp = get_xy_for_display(snap_idx, key)
                pts_pf = np.column_stack([a_pf, b_pf])
                m = pareto_mask(pts_pf)
                ax.scatter(a_disp[~m], b_disp[~m], s=18, color="0.70", alpha=0.55,
                           edgecolor="none", label="dominated")
                ax.scatter(a_disp[m],  b_disp[m],  s=44, color="#d62728",
                           edgecolor="white", lw=0.7, label="Pareto front", zorder=3)
                pf = np.column_stack([a_disp[m], b_disp[m]])
                if len(pf) > 1:
                    pf = pf[pf[:, 0].argsort()]
                    ax.plot(pf[:, 0], pf[:, 1], color="#d62728", lw=1.2, alpha=0.6, zorder=2)
                xl, yl = pair_labels[key]
                ax.set_xlabel(xl); ax.set_ylabel(yl)
                ax.set_xlim(xlim); ax.set_ylim(ylim)
                if key == "SA-Vina":
                    ax.invert_xaxis()
                elif key == "QED-SA":
                    ax.invert_yaxis()
                ax.set_title(f"{tag} bank — {key}    |PF| = {int(m.sum())}", fontsize=10)
                if row == 0 and col == 0:
                    ax.legend(loc="best", frameon=False, fontsize=8)

        # Row 2: hypervolume
        ax_hv = fig.add_subplot(gs[2, :])
        all_hv_vals = []
        for key, color, marker in [("QED-Vina", "#1f77b4", "o"),
                                    ("SA-Vina",  "#2ca02c", "s"),
                                    ("QED-SA",   "#ff7f0e", "^")]:
            vals = np.array([r[f"hv_{key.replace('-', '_')}"] for r in pareto_rows])
            ax_hv.plot(steps, vals, color=color, marker=marker, lw=2, label=key)
            all_hv_vals.append(vals)
        ax_hv.set_title("Pareto hypervolume (each metric rescaled to run-wide [0, 1])")
        ax_hv.set_xlabel("Update step")
        ax_hv.set_ylabel("Hypervolume (normalized)")
        ax_hv.set_xticks(steps)
        ax_hv.set_xlim(steps[0] - 1.5, steps[-1] + 0.5)
        hv_concat = np.concatenate(all_hv_vals)
        y_lo, y_hi = float(np.nanmin(hv_concat)), float(np.nanmax(hv_concat))
        y_span = y_hi - y_lo if y_hi > y_lo else max(abs(y_hi), 1e-3)
        pad = y_span * 0.08
        ax_hv.set_ylim(max(0.0, y_lo - pad), y_hi + pad)
        ax_hv.legend(ncol=1, frameon=False, loc="upper left")
        self._round_lines(ax_hv)

        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)


    # ============ Reaction Dynamics ============ #
    def fig_reactions(self, rxn_shares, contribution, rxn_delta, save_path):
        steps = rxn_shares.index.values
        rxns = list(rxn_shares.columns)
        colors = categorical_palette(len(rxns))

        fig = plt.figure(figsize=(13, 13), constrained_layout=False)
        fig.subplots_adjust(left=0.18, right=0.82, top=0.96, bottom=0.05,
                            hspace=0.32)
        gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.2])

        # (a) Stacked share over time
        ax1 = fig.add_subplot(gs[0])
        ax1.stackplot(steps,
                      [rxn_shares[c].values for c in rxns],
                      labels=rxns, colors=colors, alpha=0.92,
                      edgecolor="white", linewidth=0.3)
        ax1.set_title("(a) Reaction share", pad=10)
        ax1.set_ylabel("Share"); ax1.set_xlabel("Update step")
        ax1.set_ylim(0, 1); ax1.set_xticks(steps)
        ax1.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                   ncol=1, fontsize=8, frameon=False, title="Reactions")
        self._round_lines(ax1)

        # (b) Final − initial share
        ax2 = fig.add_subplot(gs[1])
        d = rxn_delta.copy()
        bar_colors = ["#2ca02c" if v > 0 else "#d62728" for v in d["delta"].values]
        y = np.arange(len(d))
        ax2.barh(y, d["delta"].values, color=bar_colors, edgecolor="white", lw=0.5)
        ax2.axvline(0, color="black", lw=0.8)
        ax2.set_yticks(y); ax2.set_yticklabels(d.index, fontsize=9)
        ax2.set_xlabel("Final share − initial share")
        ax2.set_title("(b) Reaction share change (final − initial)", pad=10)

        # (c) Top-10% enrichment
        ax3 = fig.add_subplot(gs[2])
        c = contribution.copy().sort_values("top_enrichment", ascending=True)
        n_rxn = len(c)
        y = np.arange(n_rxn)
        enrich_vals = c["top_enrichment"].values

        norm_idx = np.linspace(0, 1, n_rxn)
        cmap = plt.get_cmap("RdYlGn")
        colors_rxn = cmap(norm_idx)

        for i, val in enumerate(enrich_vals):
            ax3.hlines(y[i], 0, val, color=colors_rxn[i], lw=1.6, alpha=0.75)
        ax3.scatter(enrich_vals, y, s=110, color=colors_rxn,
                    edgecolor="black", lw=0.7, zorder=3)
        ax3.axvline(1.0, color="black", lw=1, ls="--", alpha=0.55,
                    label="enrichment = 1 (no preference)")

        delta_map = rxn_delta["delta"].to_dict() if rxn_delta is not None else {}
        for i, val in enumerate(enrich_vals):
            rname = c.index[i]
            d = float(delta_map.get(rname, 0.0))
            sign = "+" if d > 0 else ""
            d_color = "#2ca02c" if d > 0 else ("#d62728" if d < 0 else "0.5")
            ax3.text(val + 0.07, y[i], f"{val:.2f}×",
                     va="center", fontsize=9, fontweight="bold")
            ax3.text(val + 0.40, y[i], f"({sign}{d:.3f})",
                     va="center", fontsize=7.5, style="italic", color=d_color)

        ax3.set_yticks(y)
        ax3.set_yticklabels(c.index, fontsize=9)
        ax3.set_xlabel("Top-10% enrichment ratio")
        ax3.set_title("(c) Reaction quality (top-10% enrichment) and usage", pad=10)
        x_max = float(enrich_vals.max())
        ax3.set_xlim(-0.05, x_max + 1.9)
        ax3.legend(loc="lower center", frameon=False, fontsize=9)

        usage_x = x_max + 1.35
        usage_counts = c["all_count"].values if "all_count" in c.columns else None
        if usage_counts is None:
            usage_counts = (c["all_share"].values *
                            c["all_share"].values.sum()).astype(int)

        ax3.text(usage_x, n_rxn + 0.4, "n",
                 ha="center", va="center", fontsize=10, fontweight="bold",
                 color="0.25")
        for i, n in enumerate(usage_counts):
            ax3.text(usage_x, y[i], f"{int(n)}", ha="center", va="center",
                     fontsize=8.5, color="0.25",
                     bbox=dict(facecolor="white", edgecolor="0.7",
                               boxstyle="round,pad=0.25", lw=0.5))

        for ax_ext in (ax2, ax3):
            pos = ax_ext.get_position()
            ax_ext.set_position([pos.x0, pos.y0,
                                 (0.95 - pos.x0), pos.height])

        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ============ Reaction–Objective Coupling ============ #
    def fig_reaction_coupling(self, coupling_df, save_path, low_n_threshold=10):

        if coupling_df.empty:
            return

        zmat = coupling_df[["qed_z", "sa_z", "vina_z"]].values
        labels = list(coupling_df.index)
        n_vals = coupling_df["n"].values
        low_n_mask = n_vals < low_n_threshold
        n_rxn = len(labels)

        display_labels = [f"{lab} †" if low_n_mask[i] else lab
                          for i, lab in enumerate(labels)]

        fig = plt.figure(figsize=(11, max(4.5, 0.4 * n_rxn + 2.4)),
                         constrained_layout=False)
        fig.subplots_adjust(left=0.20, right=0.97, top=0.93, bottom=0.10,
                            wspace=0.25)
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.45])
        ax_h = fig.add_subplot(gs[0])
        ax_b = fig.add_subplot(gs[1])

        vmax = max(abs(np.nanmin(zmat)), abs(np.nanmax(zmat)), 0.3)

        im = ax_h.imshow(zmat, cmap="RdYlGn", aspect="auto",
                         vmin=-vmax, vmax=vmax)
        for i, is_low in enumerate(low_n_mask):
            if is_low:
                ax_h.add_patch(plt.Rectangle((-0.5, i - 0.5), 3, 1,
                                              facecolor="white", alpha=0.55,
                                              edgecolor="none", zorder=2))

        ax_h.set_xticks(range(3))
        ax_h.set_xticklabels(["QED", "SA", "−Vina"], fontsize=10)
        ax_h.set_yticks(range(n_rxn))
        ax_h.set_yticklabels(display_labels, fontsize=9)
        ax_h.set_title("(a) Per-reaction objective z-scores", fontsize=11)

        for i in range(n_rxn):
            for j in range(3):
                v = zmat[i, j]
                if np.isnan(v):
                    continue
                if low_n_mask[i]:
                    color = "0.45"
                    weight = "normal"
                else:
                    color = "white" if abs(v) > vmax * 0.6 else "black"
                    weight = "normal"
                ax_h.text(j, i, f"{v:+.2f}", ha="center", va="center",
                          fontsize=8.5, color=color, fontweight=weight, zorder=3)

        cbar = fig.colorbar(im, ax=ax_h, fraction=0.04, pad=0.02)
        cbar.set_label("z-score (run-wide reference)")

        y = np.arange(n_rxn)
        bar_colors = ["#cccccc" if low_n_mask[i] else "#666666"
                      for i in range(n_rxn)]
        bar_edges  = ["#d62728" if low_n_mask[i] else "white"
                      for i in range(n_rxn)]
        bar_widths = [1.4 if low_n_mask[i] else 0.5 for i in range(n_rxn)]
        for i in range(n_rxn):
            ax_b.barh(y[i], n_vals[i], color=bar_colors[i],
                      edgecolor=bar_edges[i], linewidth=bar_widths[i])
        ax_b.axvline(low_n_threshold, color="#d62728", lw=1, ls=":",
                     alpha=0.6, label=f"n = {low_n_threshold} threshold")
        ax_b.set_yticks(y); ax_b.set_yticklabels([])
        ax_b.set_xlabel("Total occurrences in run")
        ax_b.set_title("(b) Reaction usage", fontsize=11)
        for i, n_val in enumerate(n_vals):
            ax_b.text(n_val, i, f" {int(n_val)}", va="center", fontsize=8,
                      color="#d62728" if low_n_mask[i] else "black")
        ax_b.legend(loc="upper right", frameon=False, fontsize=8)

        shared_ylim = (-0.5, n_rxn - 0.5)
        ax_h.set_ylim(shared_ylim)
        ax_b.set_ylim(shared_ylim)
        ax_b.invert_yaxis()
        ax_h.invert_yaxis()

        if low_n_mask.any():
            fig.text(0.01, -0.01,
                     f"† Reactions with n < {low_n_threshold} are de-emphasized; "
                     "z-scores based on few samples should be interpreted cautiously.",
                     fontsize=8, style="italic", color="0.35")

        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ---------- main entry ---------- #

    def save_and_plot(self, save_dir=None, low_n_threshold=10):
        if not self.history:
            return None
        save_dir = save_dir or self.job_dir
        os.makedirs(save_dir, exist_ok=True)
        with mpl.rc_context(PLOT_STYLE):
            df, rxn_shares, contribution, rxn_delta = self._summary_frames()
            raw_delta, scaled, driver_df = self._objective_frames(df)
            pareto_rows = self._pareto_frames()
            coupling_df = self._reaction_coupling_frame()

            paths = {
                "objective": os.path.join(save_dir, "track1_objective_dynamics.png"),
                "pareto":    os.path.join(save_dir, "track2_pareto_dynamics.png"),
                "reactions": os.path.join(save_dir, "track3_reaction_dynamics.png"),
                "coupling":  os.path.join(save_dir, "track4_reaction_objective_coupling.png"),
            }
            self.fig_objective_dynamics(df, raw_delta, scaled, driver_df, paths["objective"])
            self.fig_pareto(pareto_rows, paths["pareto"])
            self.fig_reactions(rxn_shares, contribution, rxn_delta, paths["reactions"])
            self.fig_reaction_coupling(coupling_df, paths["coupling"],
                                       low_n_threshold=low_n_threshold)

            # CSV exports
            df.to_csv(os.path.join(save_dir, "summary.csv"), index=False)
            scaled.to_csv(os.path.join(save_dir, "objective_scaled.csv"), index=False)
            driver_df.to_csv(os.path.join(save_dir, "objective_driver.csv"), index=False)
            rxn_shares.to_csv(os.path.join(save_dir, "reaction_shares.csv"))
            rxn_delta.to_csv(os.path.join(save_dir, "reaction_delta.csv"))
            contribution.to_csv(os.path.join(save_dir, "reaction_contribution.csv"))
            coupling_df.to_csv(os.path.join(save_dir, "reaction_coupling.csv"))
            hv_df = pd.DataFrame([{
                "step": r["step"], "round": r["round"],
                "n_pareto_3d": r["n_pareto_3d"],
                "hv_QED_Vina_norm": r["hv_QED_Vina"],
                "hv_SA_Vina_norm":  r["hv_SA_Vina"],
                "hv_QED_SA_norm":   r["hv_QED_SA"],
                "hv_QED_Vina_raw":  r["hv_QED_Vina_raw"],
                "hv_SA_Vina_raw":   r["hv_SA_Vina_raw"],
                "hv_QED_SA_raw":    r["hv_QED_SA_raw"],
            } for r in pareto_rows])
            hv_df.to_csv(os.path.join(save_dir, "hypervolume.csv"), index=False)
            pd.DataFrame([{k: v for k, v in r.items() if k != "fronts"}
                          for r in pareto_rows]).to_csv(
                os.path.join(save_dir, "pareto.csv"), index=False)

            return paths