# SAC Analyzer GUI

A Python GUI tool for analyzing seismic data stored in SAC (Seismic Analysis Code) format.

This application reads SAC files, applies bandpass filtering, and visualizes the data in:
- Time domain (filtered signal)
- Frequency domain (FFT)
- Spectrogram

It also supports exporting the filtered signal to CSV for further analysis.

---

## 🔧 Features

- 📂 Load one or more `.sac` files  
- 🎚️ Apply Butterworth bandpass filter  
- 📊 Plot time series, frequency spectrum, and spectrogram  
- 💾 Export filtered signal as CSV  
- 🖱️ GUI built using `tkinter`, analysis with `obspy`, `scipy`, `matplotlib`

---

## 📦 Install Dependencies

```bash
pip install -r requirements.txt
```

---

## 🚀 Run the Application

```bash
python analyzer_gui.py
```

---

## 📁 Output

Filtered signals are saved as `.csv` files in the `exports/` directory.

---

## 🖼️ Screenshots

### Main Interface
![GUI](screenshots/gui_main.png)

### FFT Plot and Spectogram View
![FFT](screenshots/fft_plot_and_spectogram.png)

---

## 📝 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
