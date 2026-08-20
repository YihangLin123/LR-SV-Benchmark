import numpy as np


def summary_level_bootstrap(tp, fn, fp, n_iterations=1000, ci=95):
    n_true = tp + fn
    n_pred = tp + fp

    p_true_hit = tp / n_true if n_true > 0 else 0.0
    p_pred_correct = tp / (tp + fp) if (tp + fp) > 0 else 0.0

    boot_f1s = []
    boot_precisions = []
    boot_recalls = []

    for _ in range(n_iterations):
        tp_star_from_truth = np.random.binomial(n_true, p_true_hit)
        tp_star_from_pred = np.random.binomial(n_pred, p_pred_correct)
        tp_star = (tp_star_from_truth + tp_star_from_pred) // 2

        fn_star = n_true - tp_star
        fp_star = n_pred - tp_star
        if fp_star < 0:
            fp_star = 0

        pred_total = tp_star + fp_star
        precision_star = tp_star / pred_total if pred_total > 0 else 0.0
        recall_star = tp_star / n_true if n_true > 0 else 0.0

        if precision_star + recall_star > 0:
            f1_star = 2 * (precision_star * recall_star) / (precision_star + recall_star)
        else:
            f1_star = 0.0

        boot_f1s.append(f1_star)
        boot_precisions.append(precision_star)
        boot_recalls.append(recall_star)

    alpha = (100 - ci) / 2.0
    f1_ci = np.percentile(boot_f1s, [alpha, 100.0 - alpha])
    p_ci = np.percentile(boot_precisions, [alpha, 100.0 - alpha])
    r_ci = np.percentile(boot_recalls, [alpha, 100.0 - alpha])

    orig_p = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    orig_r = tp / n_true if n_true > 0 else 0.0
    orig_f1 = 2 * orig_p * orig_r / (orig_p + orig_r) if (orig_p + orig_r) > 0 else 0.0

    return {
        "Precision": (orig_p, p_ci),
        "Recall": (orig_r, r_ci),
        "F1": (orig_f1, f1_ci)
    }

dataset_raw = {
    "PacBio CCS 84×": [
        ("cuteSV2", 87, 69, 4306), ("Sniffles2", 112, 44, 24137), ("Debreak", 119, 37, 1604),
        ("SVIM", 118, 38, 2294), ("SVDF", 134, 22, 3693), ("SVHunter", 106, 50, 3293),
        ("NanoVar", 106, 50, 4108), ("SVision", 97, 59, 6209), ("Dysgu", 55, 101, 3887),
        ("SVision-pro", 64, 92, 394), ("Delly2", 144, 12, 107), ("nanomonsv", 140, 16, 119),
        ("SAVANA", 144, 12, 79), ("Severus", 137, 19, 50)
    ],
    "Oxford Nanopore 47×": [
        ("cuteSV2", 80, 76, 4332), ("Sniffles2", 112, 44, 24816), ("Debreak", 122, 34, 1698),
        ("SVIM", 117, 39, 1333), ("SVDF", 128, 28, 1709), ("SVHunter", 101, 55, 3920),
        ("NanoVar", 98, 58, 3053), ("SVision", 86, 70, 4162), ("Dysgu", 47, 109, 4455),
        ("SVision-pro", 48, 109, 67), ("Delly2", 140, 16, 102), ("nanomonsv", 137, 19, 49),
        ("SAVANA", 135, 21, 21), ("Severus", 148, 8, 59)
    ]
}

if __name__ == "__main__":
    np.random.seed(42)
    for dataset_name, tools_data in dataset_raw.items():
        print(f"\n================ {dataset_name} ================")
        print("Tool\tPrecision (95% CI)\tRecall (95% CI)\tF1-score (95% CI)")

        for tool, tp, fn, fp in tools_data:
            res = summary_level_bootstrap(tp, fn, fp, n_iterations=1000)
            p_val, p_ci = res["Precision"]
            r_val, r_ci = res["Recall"]
            f1_val, f1_ci = res["F1"]

            p_str = f"{p_val:.4f} [{p_ci[0]:.4f}, {p_ci[1]:.4f}]"
            r_str = f"{r_val:.4f} [{r_ci[0]:.4f}, {r_ci[1]:.4f}]"
            f1_str = f"{f1_val:.4f} [{f1_ci[0]:.4f}, {f1_ci[1]:.4f}]"

            print(f"{tool}\t{p_str}\t{r_str}\t{f1_str}")