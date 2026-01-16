# Spectral Denoising Frontend (Tkinter)

A lightweight desktop GUI (Tkinter) for running the spectral denoising pipeline on **CSV** inputs and **MSP** files.  

Please note: MSP files are preferred, as CSV files usually screw up the spectral mass/intesnity. If you really need to use CSV files, please contact me and I will give you a sample file/code to use CSV files. They require special handling beyound read in with Pandas!

Spepctral denoising removes noise ions from MS/MS spectral, based on chemical plausibility determined by structure/formula.

The original paper is published on _Nature Methods._ (https://www.nature.com/articles/s41592-025-02646-x)

---

## Features

- Desktop GUI built with Tkinter (Windows-friendly)
- Load your MSP files.
- Denoising MS/MS spectra.
- Parallel processing enabled for efficiency.

---

## Requirements

- Minimal Conda experience is required
- Python 3.12 is recommended (tested)

---

## Quick start

---

## Clone the repo

```bash
git clone <your-repo-url>
cd <your-repo-folder>
```

## Setup virtual environment
```bash
conda create -n spectral_denoising python=3.12
```

## Install dependencies
```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Finally, you can simply run the app via:
```bash
python app.py
```
**To use the app...**

### File loading

First, choose your csv file in input file tab, then click load file. This is crucial since it will not only load the file into memory but also check for file sanity.

If something is wrong, it will prompt the user to reload the file. Most .MSP files will work fine.

### Processing parameters

This is the mass error you want to use when denoising the spectra. For QE instruments, usually 0.005 Da is recommended. For ToF instruments, 0.01 is recommended.

User can also choose to use other values if desired.

### Output File

Choose where you would like to store your output, cleaned spectra file, in either MSP/CSV format. Please note  MSP files are usually recommended as CSV files require special handling.

### Process

Hit process, you are good to go! The program is usually very fast but if you have a large file (e.g. 28k spectra in VF library), it could take some time. It did not freeze!!

***Ok that's it! If you have question please reach out to me (fzkong@ucdavis.edu), or raise an issue on GitHub.***
