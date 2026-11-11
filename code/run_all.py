"""Regenerate Figs. 1-4, check them against Nature figure specifications, and print all model values."""
import subprocess, sys
for s in ["fig1.py", "fig2.py", "fig3.py", "fig4.py"]:
    subprocess.run([sys.executable, s], check=True)
for n in range(1, 5):
    subprocess.run([sys.executable, "check_nature.py", f"../figures/Fig{n}.svg"], check=True)
subprocess.run([sys.executable, "model_values.py"], check=True)
