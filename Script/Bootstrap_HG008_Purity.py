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
    "HG008 PacBio CCS 30× (100%)": [
        ("SVsision-pro", 51, 105, 111), ("Delly2", 129, 27, 45), ("nanomonsv", 131, 25, 89), ("SAVANA", 132, 24, 26), ("Severus", 136, 20, 41)
    ],
    "HG008 PacBio CCS 30× (80%)": [
        ("SVsision-pro", 51, 105, 132), ("Delly2", 131, 25, 48), ("nanomonsv", 131, 25, 84), ("SAVANA", 132, 24, 22), ("Severus", 136, 20, 38)
    ],
    "HG008 PacBio CCS 30× (60%)": [
        ("SVsision-pro", 51, 105, 169), ("Delly2", 131, 25, 48), ("nanomonsv", 130, 26, 77), ("SAVANA", 132, 24, 21), ("Severus", 136, 20, 38)
    ],
    "HG008 PacBio CCS 30× (40%)": [
        ("SVsision-pro", 51, 105, 187), ("Delly2", 131, 25, 50), ("nanomonsv", 129, 27, 67), ("SAVANA", 133, 23, 19), ("Severus", 136, 20, 31)
    ],
    "HG008 PacBio CCS 30× (20%)": [
        ("SVsision-pro", 27, 129, 143), ("Delly2", 122, 34, 25), ("nanomonsv", 121, 35, 41), ("SAVANA", 116, 40, 14), ("Severus", 129, 27, 17)
    ],
    "HG008 Oxford Nanopore 30× (100%)": [
        ("SVsision-pro", 41, 115, 22), ("Delly2", 132, 24, 64), ("nanomonsv", 134, 22, 46), ("SAVANA", 131, 25, 17), ("Severus", 144, 12, 64)
    ],
    "HG008 Oxford Nanopore 30× (80%)": [
        ("SVsision-pro", 41, 115, 23), ("Delly2", 133, 23, 61), ("nanomonsv", 134, 22, 44), ("SAVANA", 130, 26, 15), ("Severus", 143, 13, 61)
    ],
    "HG008 Oxford Nanopore 30× (60%)": [
        ("SVsision-pro", 49, 107, 18), ("Delly2", 132, 24, 57), ("nanomonsv", 133, 23, 34), ("SAVANA", 126, 30, 11), ("Severus", 141, 15, 44)
    ],
    "HG008 Oxford Nanopore 30× (40%)": [
        ("SVsision-pro", 29, 127, 17), ("Delly2", 124, 32, 32), ("nanomonsv", 127, 29, 14), ("SAVANA", 107, 49, 9), ("Severus", 139, 17, 31)
    ],
    "HG008 Oxford Nanopore 30× (20%)": [
        ("SVsision-pro", 3, 153, 16), ("Delly2", 72, 84, 14), ("nanomonsv", 93, 63, 8), ("SAVANA", 68, 88, 3), ("Severus", 108, 48, 13)
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