import pandas as pd
import matplotlib.pyplot as plt
import os

SEC_FILE = "secant_analysis.csv"
NWT_FILE = "newton_analysis.csv"

# Check if files exist
if not os.path.exists(SEC_FILE) or not os.path.exists(NWT_FILE):
    print(f"[ERROR] Required CSV files not found. Ensure {SEC_FILE} and {NWT_FILE} exist.")
    exit(1)

# Load the separate dataframes
df_sec = pd.read_csv(SEC_FILE)
df_nwt = pd.read_csv(NWT_FILE)


fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("Convergence Analysis: Secant vs Newton-Raphson", 
             fontsize=14, fontweight="bold")

# Plot 1: Successive Error (log scale)
ax = axes[0]
ax.semilogy(df_sec["iteration"], df_sec["error"], "o-", 
            color="Blue", label="Secant", linewidth=2, markersize=4)
ax.semilogy(df_nwt["iteration"], df_nwt["error"], "s-", 
            color="Magenta", label="Newton-Raphson", linewidth=2, markersize=4)

ax.set_title("Successive Error $|x_{n+1} - x_n|$")
ax.set_xlabel("Iteration")
ax.set_ylabel("Error (Log Scale)")
ax.legend()
ax.grid(True, which="both", linestyle="--", alpha=0.5)

# Plot 2: Convergence Order
ax2 = axes[1]
# We filter out the rows where order is 0.0 (usually the last iteration)
df_sec_valid = df_sec[df_sec["order"] > 0]
df_nwt_valid = df_nwt[df_nwt["order"] > 0]

ax2.plot(df_sec_valid["iteration"], df_sec_valid["order"], "o-", 
         color="Blue", label="Secant", linewidth=2, markersize=4)
ax2.plot(df_nwt_valid["iteration"], df_nwt_valid["order"], "s-", 
         color="Magenta", label="Newton-Raphson", linewidth=2, markersize=4)

# Reference lines for expected theoretical orders
ax2.axhline(y=2.0, color="green", linestyle="--", alpha=0.6, label="Order 2 (Quadratic)(Theoretical)")
ax2.axhline(y=1.618, color="orange", linestyle="--", alpha=0.6, label="Order 1.618 (Superlinear)(Theoretical)")

ax2.set_title("Estimated Convergence Order $p_n$")
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Order")
ax2.set_ylim(0, 3.5) # Keeping the scale relevant
ax2.legend(fontsize=9, loc='lower right')
ax2.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("convergence_plot_separate.png", dpi=150, bbox_inches="tight")
plt.show()
