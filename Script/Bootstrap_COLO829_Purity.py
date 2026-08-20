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
    "COLO829 PacBio CCS 30× (100%)": [
        ("SVsision-pro", 17, 51, 110), ("Delly2", 45, 23, 34), ("nanomonsv", 50, 18, 93), ("SAVANA", 47, 21, 23), ("Severus", 53, 15, 41)
    ],
    "COLO829 PacBio CCS 30× (80%)": [
        ("SVsision-pro", 17, 51, 160), ("Delly2", 47, 21, 30), ("nanomonsv", 48, 20, 58), ("SAVANA", 45, 23, 19), ("Severus", 53, 15, 41)
    ],
    "COLO829 PacBio CCS 30× (60%)": [
        ("SVsision-pro", 17, 51, 219), ("Delly2", 48, 20, 32), ("nanomonsv", 47, 21, 44), ("SAVANA", 44, 24, 19), ("Severus", 53, 15, 34)
    ],
    "COLO829 PacBio CCS 30× (40%)": [
        ("SVsision-pro", 17, 51, 210), ("Delly2", 47, 21, 26), ("nanomonsv", 42, 26, 30), ("SAVANA", 43, 25, 17), ("Severus", 53, 15, 33)
    ],
    "COLO829 PacBio CCS 30× (20%)": [
        ("SVsision-pro", 7, 61, 118), ("Delly2", 32, 36, 10), ("nanomonsv", 36, 32, 13), ("SAVANA", 30, 38, 7), ("Severus", 45, 23, 14)
    ],
    "COLO829 Oxford Nanopore 30× (100%)": [
        ("SVsision-pro", 17, 51, 40), ("Delly2", 43, 25, 67), ("nanomonsv", 46, 22, 45), ("SAVANA", 48, 20, 32), ("Severus", 50, 18, 61)
    ],
    "COLO829 Oxford Nanopore 30× (80%)": [
        ("SVsision-pro", 18, 50, 49), ("Delly2", 45, 23, 67), ("nanomonsv", 45, 23, 38), ("SAVANA", 47, 21, 26), ("Severus", 50, 18, 61)
    ],
    "COLO829 Oxford Nanopore 30× (60%)": [
        ("SVsision-pro", 17, 51, 52), ("Delly2", 44, 24, 63), ("nanomonsv", 53, 25, 30), ("SAVANA", 46, 22, 20), ("Severus", 50, 18, 53)
    ],
    "COLO829 Oxford Nanopore 30× (40%)": [
        ("SVsision-pro", 12, 56, 41), ("Delly2", 41, 27, 36), ("nanomonsv", 40, 28, 18), ("SAVANA", 39, 29, 11), ("Severus", 49, 19, 24)
    ],
    "COLO829 Oxford Nanopore 30× (20%)": [
        ("SVsision-pro", 6, 62, 27), ("Delly2", 25, 43, 12), ("nanomonsv", 25, 43, 5), ("SAVANA", 20, 48, 2), ("Severus", 35, 33, 15)
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