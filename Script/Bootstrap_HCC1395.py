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
    "PacBio CCS 62×": [
        ("cuteSV2", 631, 1187, 5327), ("Sniffles2", 516, 1302, 2577), ("Debreak", 778, 1040, 5598),
        ("SVIM", 688, 1130, 4939), ("SVDF", 920, 898, 136451), ("SVHunter", 589, 1229, 3672),
        ("NanoVar", 576, 1242, 3890), ("SVision", 548, 1270, 3813), ("Dysgu", 505, 1313, 4805),
        ("SVision-pro", 50, 1768, 2848), ("Delly2", 656, 1162, 1274), ("nanomonsv", 645, 1173, 742),
        ("SAVANA", 644, 1174, 627), ("Severus", 689, 1129, 765)
    ],
    "Oxford Nanopore 19×": [
        ("cuteSV2", 267, 1551, 8542), ("Sniffles2", 185, 1633, 7256), ("Debreak", 174, 1644, 6740),
        ("SVIM", 224, 1594, 9152), ("SVDF", 288, 1530, 47129), ("SVHunter", 218, 1600, 10596),
        ("NanoVar", 200, 1618, 6571), ("SVision", 143, 1675, 5114), ("Dysgu", 168, 1650, 8257),
        ("SVision-pro", 12, 1806, 671), ("Delly2", 1, 1817, 626), ("nanomonsv", 14, 1804, 1911),
        ("SAVANA", 1, 1817, 262), ("Severus", 30, 1788, 2971)
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