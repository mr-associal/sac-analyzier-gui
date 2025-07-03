"""
A Tkinter-based GUI application for analyzing and filtering SAC files.
Code is modularized with repeated operations turned into functions.
"""
import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from scipy.signal import butter, filtfilt, spectrogram
import os
import sys

# --- Helper Functions ---
def read_sac_file(file_path):
    """Read SAC file and return data and sampling rate."""
    st = read(file_path)
    tr = st[0]
    return tr.data, tr.stats.sampling_rate


def bandpass_filter(data, fs, lowcut, highcut, order=4):
    """
    Apply Butterworth bandpass filter based on Nyquist frequency.
    Returns original data if frequency limits are invalid.
    """
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    if not (0 < low < 1) or not (0 < high < 1) or not (low < high):
        messagebox.showwarning("Invalid Filter", f"Filtering failed.\nfs={fs}, low={lowcut}, high={highcut}")
        print(f"Warning: Invalid filter parameters! fs={fs}, low={lowcut}, high={highcut}")
        return data
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data)


def compute_fft(y_filtered, fs):
    """Compute FFT and frequency axis."""
    Y = np.fft.fft(y_filtered)
    freq = np.fft.fftfreq(len(y_filtered), d=1/fs)
    half = len(freq) // 2
    return freq[:half], np.abs(Y[:half]) / half


def compute_spectrogram(y_filtered, fs):
    """Compute spectrogram."""
    f, tt, Sxx = spectrogram(y_filtered, fs)
    return f, tt, Sxx


def get_time_axis(y, fs):
    """Create time axis."""
    T = len(y) / fs
    return np.linspace(0, T, len(y))

# --- Main Analysis Function ---
def analyze_sac_files(file_paths):
    """Analyze SAC files within filter range and plot graphs."""
    for widget in graph_frame.winfo_children():
        widget.destroy()
    fig, axs = plt.subplots(3, 1, figsize=(8, 9), dpi=100)
    try:
        fmin = float(entry_fmin.get())
        fmax = float(entry_fmax.get())
    except ValueError:
        messagebox.showerror("Error", "fmin and fmax must be numeric.")
        return
    for file_path in file_paths:
        y, fs = read_sac_file(file_path)
        t = get_time_axis(y, fs)
        y_filtered = bandpass_filter(y, fs, fmin, fmax)
        freq, Y = compute_fft(y_filtered, fs)
        axs[0].plot(t, y_filtered, label=os.path.basename(file_path))
        axs[1].plot(freq, Y, label=os.path.basename(file_path))
        f, tt, Sxx = compute_spectrogram(y_filtered, fs)
        axs[2].pcolormesh(tt, f, 10 * np.log10(Sxx), shading='gouraud')
    axs[0].set_title("Time Domain (Filtered)")
    axs[0].legend()
    axs[1].set_title("Frequency Domain (FFT)")
    axs[1].legend()
    axs[2].set_title("Spectrogram")
    fig.tight_layout()
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()

# --- Export CSV ---
def export_results():
    if not last_files:
        messagebox.showwarning("Warning", "Please select files first.")
        return
    try:
        fmin = float(entry_fmin.get())
        fmax = float(entry_fmax.get())
    except ValueError:
        messagebox.showerror("Error", "fmin and fmax must be numeric.")
        return
    os.makedirs("exports", exist_ok=True)
    for file_path in last_files:
        y, fs = read_sac_file(file_path)
        t = get_time_axis(y, fs)
        y_filtered = bandpass_filter(y, fs, fmin, fmax)
        fname = os.path.basename(file_path).replace(".sac", "_filtered.csv")
        out_path = os.path.join("exports", fname)
        np.savetxt(out_path, np.column_stack((t, y_filtered)), delimiter=",", header="Time,Amplitude")
    messagebox.showinfo("Exported", "Data saved to 'exports/' folder.")

# --- File Selection and Reset Filter ---
def browse_files():
    global last_files, last_fs, last_nyquist
    file_paths = filedialog.askopenfilenames(filetypes=[("SAC Files", "*.sac")])
    if file_paths:
        last_files = file_paths
        entry_files.delete(0, tk.END)
        entry_files.insert(0, ", ".join([os.path.basename(fp) for fp in file_paths]))
        try:
            y, fs = read_sac_file(file_paths[0])
            nyquist = fs / 2
            entry_fmin.delete(0, tk.END)
            entry_fmin.insert(0, "0.1")
            entry_fmax.delete(0, tk.END)
            entry_fmax.insert(0, f"{round(nyquist - 0.1, 2)}")
            last_fs = fs
            last_nyquist = nyquist
        except Exception as e:
            entry_fmin.delete(0, tk.END)
            entry_fmax.delete(0, tk.END)
            print(f"File reading error: {e}")
        analyze_sac_files(file_paths)


def reset_filter():
    if last_fs is None or not last_files:
        messagebox.showwarning("Warning", "No file loaded yet.")
        return
    entry_fmin.delete(0, tk.END)
    entry_fmin.insert(0, "0.1")
    entry_fmax.delete(0, tk.END)
    entry_fmax.insert(0, f"{round(last_nyquist - 0.1, 2)}")
    analyze_sac_files(last_files)

# --- GUI Setup ---
last_files = []
last_fs = None
last_nyquist = None

root = tk.Tk()
root.title("SAC Data Analyzer GUI")
root.configure(bg="#1e272e")
root.protocol("WM_DELETE_WINDOW", root.quit)

font = ("Arial", 11, "bold")
label_fg = "#ffffff"

file_frame = tk.Frame(root, bg="#1e272e")
file_frame.pack(pady=10)

entry_files = tk.Entry(file_frame, width=70)
entry_files.grid(row=0, column=0, padx=10)

tk.Button(file_frame, text="Browse SAC files", command=browse_files, font=font).grid(row=0, column=1)
tk.Button(file_frame, text="Save Data as CSV", command=export_results, font=font).grid(row=0, column=2, padx=10)
tk.Button(file_frame, text="Reset Filter", command=reset_filter, font=font).grid(row=1, column=4, padx=10)

tk.Label(file_frame, text="fmin (Hz):", bg="#1e272e", fg=label_fg, font=font).grid(row=1, column=0, pady=5, sticky='e')
entry_fmin = tk.Entry(file_frame, width=10)
entry_fmin.grid(row=1, column=1, sticky='w')

tk.Label(file_frame, text="fmax (Hz):", bg="#1e272e", fg=label_fg, font=font).grid(row=1, column=2, sticky='e')
entry_fmax = tk.Entry(file_frame, width=10)
entry_fmax.grid(row=1, column=3, sticky='w')

graph_frame = tk.Frame(root, bg="#1e272e")
graph_frame.pack()

root.mainloop()
sys.exit()