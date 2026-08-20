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
    "Oxford Nanopore 49×": [
        ("cuteSV2", 41, 27, 6536), ("Sniffles2", 37, 31, 1135), ("Debreak", 43, 25, 1433),
        ("SVIM", 41, 27, 1158), ("SVDF", 49, 19, 4842), ("SVHunter", 24, 44, 3761),
        ("NanoVar", 35, 33, 2550), ("SVision", 29, 39, 2802), ("Dysgu", 14, 54, 4114),
        ("SVision-pro", 18, 50, 66), ("Delly2", 46, 22, 90), ("nanomonsv", 46, 22, 48),
        ("SAVANA", 50, 18, 38), ("Severus", 49, 19, 64)
    ],
    "Oxford Nanopore 40×": [
        ("cuteSV2", 21, 47, 2230), ("Sniffles2", 37, 31, 1205), ("Debreak", 30, 38, 982),
        ("SVIM", 35, 33, 1127), ("SVDF", 46, 22, 2488), ("SVHunter", 23, 45, 3777),
        ("NanoVar", 34, 34, 2889), ("SVision", 30, 38, 2949), ("Dysgu", 13, 55, 3831),
        ("SVision-pro", 18, 50, 62), ("Delly2", 45, 23, 85), ("nanomonsv", 48, 23, 49),
        ("SAVANA", 50, 18, 33), ("Severus", 50, 18, 68)
    ],
    "Oxford Nanopore 30×": [
        ("cuteSV2", 21, 47, 2182), ("Sniffles2", 38, 30, 1311), ("Debreak", 37, 31, 1151),
        ("SVIM", 34, 34, 1406), ("SVDF", 41, 27, 3151), ("SVHunter", 26, 42, 4104),
        ("NanoVar", 32, 36, 3596), ("SVision", 27, 41, 3740), ("Dysgu", 13, 55, 3440),
        ("SVision-pro", 17, 51, 40), ("Delly2", 43, 25, 67), ("nanomonsv", 46, 22, 45),
        ("SAVANA", 48, 20, 32), ("Severus", 50, 18, 61)
    ],
    "Oxford Nanopore 20×": [
        ("cuteSV2", 19, 49, 2795), ("Sniffles2", 32, 36, 1399), ("Debreak", 33, 35, 1509),
        ("SVIM", 34, 34, 2154), ("SVDF", 27, 41, 4893), ("SVHunter", 22, 46, 4560),
        ("NanoVar", 22, 46, 5186), ("SVision", 15, 53, 6048), ("Dysgu", 12, 56, 2818),
        ("SVision-pro", 12, 56, 33), ("Delly2", 39, 29, 51), ("nanomonsv", 46, 22, 23),
        ("SAVANA", 41, 27, 26), ("Severus", 48, 20, 56)
    ],
    "Oxford Nanopore 10×": [
        ("cuteSV2", 19, 49, 7996), ("Sniffles2", 24, 44, 3149), ("Debreak", 19, 49, 4181),
        ("SVIM", 21, 47, 5055), ("SVDF", 24, 44, 4935), ("SVHunter", 16, 49, 5867),
        ("NanoVar", 24, 44, 5492), ("SVision", 16, 52, 7170), ("Dysgu", 11, 57, 2833),
        ("SVision-pro", 21, 47, 84), ("Delly2", 14, 54, 46), ("nanomonsv", 32, 36, 72),
        ("SAVANA", 24, 44, 14), ("Severus", 37, 31, 109)
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