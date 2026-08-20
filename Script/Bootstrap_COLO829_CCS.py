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
        if fp_star < 0: fp_star = 0

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
    "PacBio CCS 66×": [
        ("cuteSV2", 49, 19, 12602), ("Sniffles2", 37, 31, 849), ("Debreak", 45, 23, 1986),
        ("SVIM", 52, 16, 1416), ("SVDF", 60, 8, 37563), ("SVHunter", 29, 39, 2680),
        ("NanoVar", 35, 33, 2056), ("SVision", 32, 36, 2262), ("Dysgu", 13, 55, 2706),
        ("SVision-pro", 24, 44, 217), ("Delly2", 51, 17, 61), ("nanomonsv", 50, 18, 40),
        ("SAVANA", 53, 15, 23), ("Severus", 55, 13, 40)
    ],
    "PacBio CCS 50×": [
        ("cuteSV2", 44, 24, 4210), ("Sniffles2", 42, 26, 2455), ("Debreak", 43, 25, 1513),
        ("SVIM", 50, 18, 1543), ("SVDF", 52, 16, 2088), ("SVHunter", 27, 41, 2765),
        ("NanoVar", 34, 34, 2311), ("SVision", 33, 35, 1624), ("Dysgu", 13, 55, 2753),
        ("SVision-pro", 22, 46, 181), ("Delly2", 50, 18, 42), ("nanomonsv", 50, 18, 53),
        ("SAVANA", 53, 15, 19), ("Severus", 54, 14, 42)
    ],
    "PacBio CCS 40×": [
        ("cuteSV2", 38, 30, 4238), ("Sniffles2", 42, 26, 2592), ("Debreak", 39, 29, 1824),
        ("SVIM", 50, 18, 1724), ("SVDF", 49, 19, 2225), ("SVHunter", 27, 41, 2822),
        ("NanoVar", 33, 35, 2402), ("SVision", 31, 37, 4406), ("Dysgu", 14, 54, 2601),
        ("SVision-pro", 20, 48, 148), ("Delly2", 50, 18, 39), ("nanomonsv", 50, 18, 84),
        ("SAVANA", 53, 15, 19), ("Severus", 53, 15, 45)
    ],
    "PacBio CCS 30×": [
        ("cuteSV2", 34, 34, 4764), ("Sniffles2", 41, 27, 2877), ("Debreak", 35, 33, 2690),
        ("SVIM", 45, 23, 2048), ("SVDF", 43, 25, 2670), ("SVHunter", 24, 44, 3079),
        ("NanoVar", 30, 38, 2735), ("SVision", 26, 42, 2999), ("Dysgu", 13, 55, 2353),
        ("SVision-pro", 17, 51, 110), ("Delly2", 45, 23, 34), ("nanomonsv", 50, 18, 93),
        ("SAVANA", 47, 21, 23), ("Severus", 53, 15, 41)
    ],
    "PacBio CCS 20×": [
        ("cuteSV2", 23, 45, 3181), ("Sniffles2", 31, 37, 1459), ("Debreak", 32, 36, 1924),
        ("SVIM", 36, 32, 2129), ("SVDF", 39, 29, 2447), ("SVHunter", 23, 45, 3288),
        ("NanoVar", 32, 36, 3162), ("SVision", 29, 39, 3268), ("Dysgu", 12, 56, 2252),
        ("SVision-pro", 20, 48, 369), ("Delly2", 38, 30, 38), ("nanomonsv", 46, 22, 103),
        ("SAVANA", 44, 24, 19), ("Severus", 49, 19, 46)
    ],
    "PacBio CCS 10×": [
        ("cuteSV2", 16, 52, 5475), ("Sniffles2", 23, 45, 2832), ("Debreak", 18, 50, 3490),
        ("SVIM", 25, 43, 3634), ("SVDF", 26, 42, 3898), ("SVHunter", 23, 45, 3898),
        ("NanoVar", 18, 50, 4042), ("SVision", 15, 53, 4494), ("Dysgu", 10, 58, 2414),
        ("SVision-pro", 12, 56, 277), ("Delly2", 13, 55, 48), ("nanomonsv", 37, 31, 144),
        ("SAVANA", 24, 44, 14), ("Severus", 39, 29, 97)
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