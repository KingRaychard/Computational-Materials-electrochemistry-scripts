import numpy as np
import matplotlib.pyplot as plt
# -----------------------------
# TO GENERATE OUTPUT ENTER THE FILENAME (.csv) WHEN RUNNING PROGRAM
# ex. /Users/sohamkulkarni/Li2.csv
# -----------------------------
filename = input("Enter path to .csv file (MUST BE FULL PATH): ")
data = []
with open(filename, "r", errors="ignore") as f:
    reading = False
    for line in f:
        line = line.strip()
        if line == "Time/sec, Current/A":
            reading = True
            continue
        if reading and line:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 2:
                try:
                    t = float(parts[0])
                    i = float(parts[1])
                    data.append((t, i))
                except ValueError:
                    pass
if not data:
    raise ValueError("No data were parsed from the file.")
times = np.array([d[0] for d in data], dtype=float)
currents = np.array([d[1] for d in data], dtype=float)

# Clean bad / repeated time rows
valid = np.r_[True, np.diff(times) > 0]
times = times[valid]
currents = currents[valid]
if len(times) < 2:
    raise ValueError("Not enough valid increasing time points after cleaning.")
# Safety check
if np.any(np.diff(times) <= 0):
    raise ValueError("Time values are not strictly increasing.")
# -----------------------------
# Find pulse boundaries
# The file switches current direction roughly every 120 s.
# Large jumps in current mark the pulse transitions.
# -----------------------------
jump_threshold = 1e-3 # A; works well for this file
jump_idx = np.where(np.abs(np.diff(currents)) > jump_threshold)[0]
# Segment boundaries
bounds = np.r_[0, jump_idx + 1, len(times)]
pulse_charges_mC = []
pulse_mean_currents_mA = []
# -----------------------------
# Plot 1: cumulative charge within each pulse
# This lets you compare all pulses directly.
# -----------------------------
for n, (a, b) in enumerate(zip(bounds[:-1], bounds[1:]), start=1):
    t_seg = times[a:b]
    i_seg = currents[a:b]
    if len(t_seg) < 2:
        continue
    # Reset time inside each pulse
    t_rel = t_seg - t_seg[0]
    # Trapezoidal cumulative integration
    dt = np.diff(t_seg)
    avg_i = 0.5 * (i_seg[:-1] + i_seg[1:])
    dq_mC = avg_i * dt * 1000.0 # A*s -> C -> mC
    q_cum_mC = np.r_[0.0, np.cumsum(dq_mC)]
    # Final pulse charge
    pulse_charges_mC.append(q_cum_mC[-1])
    pulse_mean_currents_mA.append(np.mean(i_seg) * 1000.0)

# -----------------------------
# Plot 2: charge delivered per pulse
# If the sets are behaving consistently, this is the plot
# that should be near-horizontal for each type of pulse.
# -----------------------------
pulse_charges_mC = np.array(pulse_charges_mC)
pulse_mean_currents_mA = np.array(pulse_mean_currents_mA)
pulse_numbers = np.arange(1, len(pulse_charges_mC) + 1)
# Separate odd/even pulses (charge/discharge)
odd = pulse_numbers % 2 == 1
even = ~odd
plt.figure(figsize=(10, 5))
plt.plot(pulse_numbers[odd], pulse_charges_mC[odd], "o-", label="Odd pulses")
plt.plot(pulse_numbers[even], pulse_charges_mC[even], "o-", label="Even pulses")
plt.xlabel("Pulse number")
plt.ylabel("Charge per pulse (mC)")
plt.title("Charge per Pulse")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# Summary
print(f"Total duration: {times[-1] / 3600:.2f} hr")
print(f"Number of pulses found: {len(pulse_charges_mC)}")
print(f"Mean odd-pulse charge:  {np.mean(pulse_charges_mC[odd]):.2f} mC")
print(f"Mean even-pulse charge: {np.mean(pulse_charges_mC[even]):.2f} mC")
print(f"Mean odd-pulse current:  {np.mean(pulse_mean_currents_mA[odd]):.3f} mA")
print(f"Mean even-pulse current: {np.mean(pulse_mean_currents_mA[even]):.3f} mA")

