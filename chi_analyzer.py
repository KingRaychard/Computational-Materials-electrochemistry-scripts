#!/usr/bin/python3

# it should be painfully obvious that this done in collaboration with Dr. ChatGPT

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import sys


print("Thor is here\n") #Of course this can only be Ray's doing

# =========================================================
# general utilities
# =========================================================

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyze CH Instruments electrochemistry data."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input CH Instruments CSV file."
    )
    parser.add_argument(
        "--technique",
	required=True,
        choices=[
            "cyclicvoltammetry",
        ],
	help=(
	    "implemented options are",
	    "cyclicvoltammetry"
        )
    )
    parser.add_argument(
        "--show_cycles",
        nargs="*",
        metavar="CYCLE",
        help=(
            "Plot CV cycles. With no cycle numbers, plot all cycles. "
            "Specify individual cycles or ranges, e.g. "
            "--show_cycles 1 3 5-8"
        )
    )
    parser.add_argument(
        "--show_combined_cycles",
        nargs="*",
        metavar="CYCLE",
        help=(
            "Plot selected CV cycles together. With no cycle numbers, "
            "plot all cycles. Individual cycles and ranges can be used, "
            "e.g. --combined_cycles 1 3 5-8"
        )
    )
    return parser.parse_args()

# =========================================================
# cyclic voltammetry
# =========================================================

def read_cyclicvoltammetry(filename):
    """
    Read cyclic voltammetry data from a CH Instruments CSV file.
    CV data starts at the line:
        Potential/V, Current/A
    Returns
    -------
    pandas.DataFrame
        DataFrame containing CV potential and current data.
    """
    cv_header = "Potential/V, Current/A"
    header_line = None
    try:
        with open(filename, "r", encoding="utf-8-sig") as f:
            for line_number, line in enumerate(f):
                if line.strip() == cv_header:
                    header_line = line_number
                    break
    except FileNotFoundError:
        print(f"Error: file not found: {filename}", file=sys.stderr)
        sys.exit(1)
    if header_line is None:
        print(
            f'Error: CV data header "{cv_header}" not found in {filename}',
            file=sys.stderr
        )
        sys.exit(1)
    try:
        data = pd.read_csv(
            filename,
            skiprows=header_line,
            skipinitialspace=True
        )
        data = data.rename(columns={
            "Potential/V": "Potential (V)",
            "Current/A": "Current (A)"
        })
    except Exception as exc:
        print(f"Error reading CV data: {exc}", file=sys.stderr)
        sys.exit(1)
    return data

def split_cv_cycles(data):
    """
    Split cyclic voltammetry data into individual cycles.
    Parameters
    ----------
    data : pandas.DataFrame
        CV data containing a "Potential (V)" column.
    Returns
    -------
    list of pandas.DataFrame
        Each element is one complete CV cycle.
    """
    potential = data["Potential (V)"].to_numpy()
    # Determine scan direction from change in potential
    dp = np.diff(potential)
    direction = np.sign(dp)
    # Deal with consecutive points having exactly the same potential.
    # Carry the previous scan direction through these points.
    for i in range(1, len(direction)):
        if direction[i] == 0:
            direction[i] = direction[i - 1]
    # Find locations where scan direction changes
    turning_points = np.where(direction[1:] != direction[:-1])[0] + 1
    cycles = []
    # A complete CV cycle contains two scan-direction changes
    start = 0
    for i in range(1, len(turning_points), 2):
        end = turning_points[i] + 1
        cycle = data.iloc[start:end].copy()
        cycle.reset_index(drop=True, inplace=True)
        cycles.append(cycle)
        # The ending point of one cycle is also the starting point
        # of the next cycle
        start = turning_points[i]
    # Add whatever remains as an incomplete final cycle
    if start < len(data) - 1:
        cycle = data.iloc[start:].copy()
        cycle.reset_index(drop=True, inplace=True)
        cycles.append(cycle)
    return cycles

def split_cv_segments(data):
    """
    Split one CV cycle into monotonic potential segments.
    """
    voltage = data["Potential (V)"].to_numpy()
    dp = np.diff(voltage)
    direction = np.sign(dp)
    for i in range(1, len(direction)):
        if direction[i] == 0:
            direction[i] = direction[i - 1]
    turning_points = np.where(direction[1:] != direction[:-1])[0] + 1
    boundaries = [0] + turning_points.tolist() + [len(data) - 1]
    segments = []
    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1] + 1
        segment = data.iloc[start:end].copy()
        segment.reset_index(drop=True, inplace=True)
        segments.append(segment)
    return segments

def calculate_relative_prominence(voltage, current, peak_index,
                                  absolute_prominence,
                                  voltage_window=0.25):
    """
    Calculate peak prominence relative to the local current range.
    voltage_window is the half-width of the desired local region.
    For example, voltage_window=0.25 uses a 0.50 V region.
    Near the edge of the scan, the window is shifted inward so that
    its total width remains approximately constant when possible.
    """
    peak_voltage = voltage[peak_index]
    voltage_min = np.min(voltage)
    voltage_max = np.max(voltage)
    full_width = 2 * voltage_window
    lower_bound = peak_voltage - voltage_window
    upper_bound = peak_voltage + voltage_window
    if lower_bound < voltage_min:
        lower_bound = voltage_min
        upper_bound = min(voltage_min + full_width, voltage_max)
    if upper_bound > voltage_max:
        upper_bound = voltage_max
        lower_bound = max(voltage_max - full_width, voltage_min)
    mask = (voltage >= lower_bound) & (voltage <= upper_bound)
    local_current = current[mask]

    if len(local_current) == 0:
        return 0.0
    local_range = np.max(local_current) - np.min(local_current)
    if local_range == 0:
        return 0.0
    return absolute_prominence / local_range


def separate_close_peaks(peaks, voltage_threshold=0.05):
    """
    Group peaks that occur within voltage_threshold of one another.
    Within each group, the peak with the largest absolute prominence
    is retained as the real peak. Other peaks are labeled as noise.
    peaks is a list of dictionaries.
    """
    if len(peaks) < 2:
        return peaks, []
    peaks = sorted(peaks, key=lambda x: x["potential"])
    groups = [[peaks[0]]]
    for peak in peaks[1:]:
        if abs(peak["potential"] - groups[-1][-1]["potential"]) < voltage_threshold:
            groups[-1].append(peak)
        else:
            groups.append([peak])
    normal_peaks = []
    noise_peaks = []
    for group in groups:
        if len(group) == 1:
            normal_peaks.append(group[0])
            continue
        best_peak = max(group, key=lambda x: x["absolute_prominence"])
        normal_peaks.append(best_peak)
        for peak in group:
            if peak is not best_peak:
                noise_peaks.append(peak)
    return normal_peaks, noise_peaks


def find_peaks_in_segment(voltage, current, peak_type,
                          voltage_window=0.25,
                          voltage_threshold=0.05,
                          min_absolute_prominence=0.0,
                          min_relative_prominence=0.0):
    """
    Find oxidation or reduction peaks within one monotonic CV segment.
    Parameters
    ----------
    peak_type : str
        "maximum" searches for positive-current maxima.
        "minimum" searches for negative-current minima.
    min_absolute_prominence : float
        Minimum absolute prominence in A.
    min_relative_prominence : float
        Minimum relative prominence expressed as a fraction.
        For example, 0.10 means 10%.
    Returns
    -------
    normal_peaks : list
    noise_peaks : list
    """

    if peak_type == "maximum":
        search_current = current
    elif peak_type == "minimum":
        search_current = -current
    else:
        raise ValueError('peak_type must be "maximum" or "minimum"')

    indices, properties = find_peaks(
        search_current,
        prominence=min_absolute_prominence
    )
    peaks = []
    for i, peak_index in enumerate(indices):
        absolute_prominence = properties["prominences"][i]

        relative_prominence = calculate_relative_prominence(
            voltage,
            search_current,
            peak_index,
            absolute_prominence,
            voltage_window
        )
        if relative_prominence < min_relative_prominence:
            continue
        peaks.append({
            "index": peak_index,
            "potential": voltage[peak_index],
            "current": current[peak_index],
            "absolute_prominence": absolute_prominence,
            "relative_prominence": relative_prominence
        })

    normal_peaks, noise_peaks = separate_close_peaks(
        peaks,
        voltage_threshold
    )
    return normal_peaks, noise_peaks


def find_cv_peaks(data, voltage_window=0.25,
                  voltage_threshold=0.05,
                  min_absolute_prominence=0.0,
                  min_relative_prominence=0.0,
                  oxidation_sign=1):
    """
    Find oxidation and reduction peaks in one CV cycle.
    oxidation_sign:
         1 = oxidation current is positive
        -1 = oxidation current is negative
    If no oxidation or reduction peak is detected, the corresponding
    current extreme and potential are returned as a fallback.
    """

    segments = split_cv_segments(data)
    oxidation_peaks = []
    reduction_peaks = []
    oxidation_noise = []
    reduction_noise = []

    for segment_number, segment in enumerate(segments, start=1):
        voltage = segment["Potential (V)"].to_numpy()
        current = segment["Current (A)"].to_numpy()

        # Positive-going scan: only search for oxidation peaks
        if voltage[-1] > voltage[0]:
            peak_type = "maximum" if oxidation_sign == 1 else "minimum"

            peaks, noise = find_peaks_in_segment(
                voltage,
                current,
                peak_type,
                voltage_window,
                voltage_threshold,
                min_absolute_prominence,
                min_relative_prominence
            )

            for peak in peaks:
                peak["segment"] = segment_number
            for peak in noise:
                peak["segment"] = segment_number

            oxidation_peaks.extend(peaks)
            oxidation_noise.extend(noise)

        # Negative-going scan: only search for reduction peaks
        elif voltage[-1] < voltage[0]:
            peak_type = "minimum" if oxidation_sign == 1 else "maximum"

            peaks, noise = find_peaks_in_segment(
                voltage,
                current,
                peak_type,
                voltage_window,
                voltage_threshold,
                min_absolute_prominence,
                min_relative_prominence
            )

            for peak in peaks:
                peak["segment"] = segment_number
            for peak in noise:
                peak["segment"] = segment_number

            reduction_peaks.extend(peaks)
            reduction_noise.extend(noise)
    return {
        "oxidation_peaks": oxidation_peaks,
        "reduction_peaks": reduction_peaks,
        "oxidation_noise": oxidation_noise,
        "reduction_noise": reduction_noise
    }

def analyze_cv(data,args):
    """
    Analyze cyclic voltammetry data.
    """
    #Split the data into cycles.
    cycles = split_cv_cycles(data)
    #plot cycles if asked to
    if args.show_cycles is not None:
        plot_cv_cycles(cycles, args.show_cycles)
    if args.show_combined_cycles is not None:
        plot_combined_cv_cycles(cycles, args.show_combined_cycles)

    print("Number of cycles: ", len(cycles))
    all_peaks = []
    for cycle in cycles:
        cycle_peaks = find_cv_peaks(
            cycle,
            voltage_window=0.25,
            voltage_threshold=0.05,
            min_absolute_prominence=0.0,
            min_relative_prominence=0.10,
            oxidation_sign=-1 # for CHI when oxidation (positive potential, 
                              # anodic current, is negative)
        )
        all_peaks.append(cycle_peaks)
    print("debug printing all oxidation peaks")
    for i in all_peaks:
        print(i["oxidation_peaks"])
    print("debug printing all oxidation noise")
    for i in all_peaks:
        print(i["oxidation_noise"])

def parse_cycle_selection(selection, number_of_cycles):
    """
    Convert cycle selections such as ["1", "3", "5-8"] into
    a sorted list of cycle numbers.
    Cycle numbering starts at 1.
    """
    if selection is None:
        return None
    if len(selection) == 0:
        return list(range(1, number_of_cycles + 1))
    selected = set()
    for item in selection:
        if "-" in item:
            try:
                start, end = map(int, item.split("-", 1))
            except ValueError:
                raise ValueError(f"Invalid cycle range: {item}")
            if start > end:
                raise ValueError(f"Invalid cycle range: {item}")
            selected.update(range(start, end + 1))
        else:
            try:
                selected.add(int(item))
            except ValueError:
                raise ValueError(f"Invalid cycle number: {item}")
    for cycle_number in selected:
        if cycle_number < 1 or cycle_number > number_of_cycles:
            raise ValueError(
                f"Cycle {cycle_number} does not exist. "
                f"Available cycles: 1-{number_of_cycles}"
            )
    return sorted(selected)

def plot_cv_cycles(cycles, selection=None):
    """
    Plot selected CV cycles in separate figures.

    Parameters
    ----------
    cycles : list
        List of CV cycle DataFrames.

    selection : list or None
        Cycle selections from --show_cycles.
        [] means show all cycles.
    """
    cycle_numbers = parse_cycle_selection(selection, len(cycles))

    if cycle_numbers is None:
        return

    for cycle_number in cycle_numbers:
        cycle = cycles[cycle_number - 1]
        voltage = cycle["Potential (V)"].to_numpy()
        current = cycle["Current (A)"].to_numpy()
        plt.figure()
        plt.plot(voltage, current)
        plt.xlabel("Potential (V)")
        plt.ylabel("Current (A)")
        plt.title(f"CV Cycle {cycle_number}")
        plt.grid(True)
        plt.gca().invert_xaxis()
    plt.show()


def plot_combined_cv_cycles(cycles, selection=None):
    """
    Plot selected CV cycles together using rainbow colors.
    The first selected cycle is red and the last selected cycle
    is purple.
    """
    cycle_numbers = parse_cycle_selection(selection, len(cycles))
    if cycle_numbers is None:
        return
    colors = plt.cm.Greys(np.linspace(0.3, 0.9, len(cycle_numbers)))
    plt.figure()
    for cycle_number, color in zip(cycle_numbers, colors):
        cycle = cycles[cycle_number - 1]
        voltage = cycle["Potential (V)"].to_numpy()
        current = cycle["Current (A)"].to_numpy()
        plt.plot(
            voltage,
            current,
            color=color,
            label=f"Cycle {cycle_number}"
        )
    plt.xlabel("Potential (V)")
    plt.ylabel("Current (A)")
    plt.title("CV Cycles")
    plt.gca().invert_xaxis()
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.gca().invert_xaxis()

# =========================================================
# main
# =========================================================

def main():
    args = parse_arguments()
    if args.technique=="cyclicvoltammetry":
        print(f"Loaded as cyclic voltammetry data: {args.input}")
        data = read_cyclicvoltammetry(args.input)
        analyze_cv(data,args)

    # For now:


if __name__ == "__main__":
    main()
