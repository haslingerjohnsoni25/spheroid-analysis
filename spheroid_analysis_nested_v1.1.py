import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
from termcolor import colored
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
import statsmodels.formula.api as smf
from sklearn.utils import resample
from skimage.filters.rank import entropy
from skimage.morphology import disk

# =============================================================================
# Calibration Constants
# =============================================================================
um_per_px = 200.0 / 298.412
area_conversion = um_per_px ** 2
print(colored("🔬 Calibration set: {:.4f} µm/px and {:.4f} µm²/px²".format(um_per_px, area_conversion), "cyan"))

# =============================================================================
# Function: Analyze Image for Spheroid Fragmentation and Complexity
# =============================================================================
def analyze_image_entropy(image_path):
    print(colored(f"🖼️  Processing: {image_path}", "yellow"))
    img_color = cv2.imread(image_path)
    if img_color is None:
        print(colored(f"⚠️  Warning: Could not load image: {image_path}", "red"))
        return None

    img = cv2.cvtColor(img_color, cv2.COLOR_BGR2GRAY)
    equalized = cv2.equalizeHist(img)
    blurred = cv2.GaussianBlur(equalized, (15, 15), 0)
    adaptive_thresh = cv2.adaptiveThreshold(blurred, 255,
                                            cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                            cv2.THRESH_BINARY_INV, 11, 2)

    # Morphological filtering
    kernel = np.ones((5, 5), np.uint8)
    cleaned = cv2.morphologyEx(adaptive_thresh, cv2.MORPH_OPEN, kernel, iterations=2)

    # Distance transform and watershed segmentation
    dist_transform = cv2.distanceTransform(cleaned, cv2.DIST_L2, 5)
    _, sure_fg = cv2.threshold(dist_transform, 0.3 * dist_transform.max(), 255, 0)
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(cleaned, sure_fg)
    _, markers = cv2.connectedComponents(sure_fg)
    markers = markers + 1
    markers[unknown == 255] = 0
    markers = cv2.watershed(cv2.cvtColor(img_color, cv2.COLOR_BGR2RGB), markers)

    # Overlay boundaries
    overlay = img_color.copy()
    overlay[markers == -1] = [0, 255, 0]  # green boundary

    # Count fragments
    num_fragments = np.max(markers) - 1  # subtract background

    # Entropy map inside detected spheroid mask
    mask = np.uint8(sure_fg)
    entropy_map = entropy(img, disk(5))
    masked_entropy = entropy_map[mask == 255]
    entropy_value = float(np.mean(masked_entropy)) if masked_entropy.size > 0 else 0

    # Contour detection for holes
    contours, _ = cv2.findContours(mask, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    num_holes = 0
    for i in range(len(contours)):
        if _[0][i][3] != -1:  # hole has a parent contour
            num_holes += 1

    # Total area
    total_area_um2 = np.sum(mask == 255) * area_conversion

    # Annotate visual
    cv2.putText(overlay, f"Fragments: {num_fragments}", (10, 25), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 2)
    cv2.putText(overlay, f"Entropy: {entropy_value:.2f}", (10, 50), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 2)
    cv2.putText(overlay, f"Holes: {num_holes}", (10, 75), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 2)

    # Save output
    out_path = os.path.splitext(image_path)[0] + "_frag_detected.png"
    cv2.imwrite(out_path, overlay)
    print(colored(f"📸 Overlay saved: {out_path}", "green"))

    return total_area_um2, num_fragments, entropy_value, num_holes

# =============================================================================
# Traverse Folders and Collect Data
# =============================================================================
def collect_entropy_from_dirs(parent_dir, group_label):
    records = []
    for well in os.listdir(parent_dir):
        well_path = os.path.join(parent_dir, well)
        if not os.path.isdir(well_path):
            continue
        print(colored(f"🔍 Scanning well: {well}", "blue"))
        for fname in os.listdir(well_path):
            if fname.lower().endswith(('.jpg', '.jpeg', '.png')):
                session = os.path.splitext(fname)[0]
                image_path = os.path.join(well_path, fname)
                result = analyze_image_entropy(image_path)
                if result:
                    area, fragments, entropy_val, holes = result
                    records.append({
                        "Image": well,
                        "Session": session,
                        "Group": group_label,
                        "Area_um2": area,
                        "Fragments": fragments,
                        "Entropy": entropy_val,
                        "Holes": holes
                    })
    return pd.DataFrame(records)

# =============================================================================
# Run Analysis
# =============================================================================
print(colored("📊 Starting advanced fragmentation analysis for all wells...", "magenta"))
control_dir = "/media/aguilarr@campus.wra.net/armoratd/prca_vaccine/prev_analysis/prevention_control"
experimental_dir = "/media/aguilarr@campus.wra.net/armoratd/prca_vaccine/prev_analysis/prevention_experimental"
control_df = collect_entropy_from_dirs(control_dir, "Control")
experimental_df = collect_entropy_from_dirs(experimental_dir, "Experimental")
combined_df = pd.concat([control_df, experimental_df], ignore_index=True)

# Format and save
combined_df["Session"] = pd.to_datetime(combined_df["Session"], errors='coerce')
combined_df.dropna(subset=["Session"], inplace=True)
combined_df.sort_values(by=["Group", "Image", "Session"], inplace=True)

excel_file = "fragmentation_results.xlsx"
combined_df.to_excel(excel_file, index=False)
print(colored(f"✅ Results saved to: {excel_file}", "green"))

# =============================================================================
# Plot Entropy by Group and Time
# =============================================================================
print(colored("📈 Generating plots...", "cyan"))
plt.figure(figsize=(10, 6))
sns.lineplot(data=combined_df, x="Session", y="Entropy", hue="Group", marker="o", ci="sd")
plt.title("Cytotoxicity (Entropy) Over Time")
plt.ylabel("Entropy (Local Texture Variation)")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("entropy_fragments_plot.png")
print(colored("📊 Entropy trend plot saved: entropy_fragments_plot.png", "cyan"))
