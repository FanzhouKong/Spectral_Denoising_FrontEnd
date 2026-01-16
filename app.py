import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import pandas as pd
import threading
from rdkit import Chem
from spectral_denoising import file_io, spectral_denoising
from pipeline import denoising_pipeline

class SpectralDenoisingGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Spectral Denoising Application")
        self.root.geometry("900x700")
        
        # Variables
        self.input_file_path = tk.StringVar()
        self.output_file_path = tk.StringVar()
        self.file_type = tk.StringVar(value="msp")
        self.mass_tolerance = tk.DoubleVar(value=0.005)
        self.mass_mode = tk.StringVar(value="orbitrap")  # orbitrap/tof/custom
        self.loaded_data = None
        
        self.create_widgets()
    
    def standardize_columns(self, df):
        """Standardize column names to a consistent format (case-insensitive mapping)"""
        # Define column mappings
        column_mapping = {
            # Peaks/spectra variations
            'msms': 'peaks',
            'ms/ms': 'peaks',
            'peak': 'peaks',
            'peaks': 'peaks',
            'spectra': 'peaks',
            'spectrum': 'peaks',
            # SMILES variations
            'smile': 'smiles',
            'smiles': 'smiles',
            # Adduct variations
            'adduct': 'adducts',
            'adducts': 'adducts',
        }
        
        # Create a mapping of current column names to standardized names
        new_columns = {}
        for col in df.columns:
            col_lower = col.lower()
            if col_lower in column_mapping:
                new_columns[col] = column_mapping[col_lower]
                self.log_message(f"Mapped column '{col}' -> '{column_mapping[col_lower]}'")
        
        # Rename columns
        if new_columns:
            df = df.rename(columns=new_columns)
            self.log_message(f"Column standardization complete")
        else:
            self.log_message("No column names needed standardization")
        
        return df
    
    def extract_smiles_from_comments(self, comments_str):
        """Extract SMILES string from comments field"""
        if pd.isna(comments_str) or not isinstance(comments_str, str):
            return None
        
        # Look for 'computed SMILES=' pattern
        import re
        match = re.search(r'computed SMILES=([^"]+)', comments_str)
        if match:
            return match.group(1).strip()
        
        # Also try 'SMILES=' pattern (without 'computed')
        match = re.search(r'SMILES=([^"]+)', comments_str)
        if match:
            return match.group(1).strip()
        
        return None
    
    def fill_missing_smiles(self, df):
        """Fill missing SMILES column from comments if available"""
        if 'smiles' not in df.columns and 'comments' in df.columns:
            self.log_message("SMILES column not found, attempting to extract from comments...")
            df['smiles'] = df['comments'].apply(self.extract_smiles_from_comments)
            
            # Count how many were successfully extracted
            valid_smiles = df['smiles'].notna().sum()
            total_rows = len(df)
            self.log_message(f"Extracted SMILES for {valid_smiles}/{total_rows} entries from comments")
            
            if valid_smiles == 0:
                self.log_message("Warning: No SMILES found in comments")
        
        return df
    
    def validate_smiles(self, df):
        """Validate all SMILES in the dataframe using RDKit"""
        if 'smiles' not in df.columns:
            return df
        
        self.log_message("Validating SMILES structures...")
        
        invalid_indices = []
        invalid_smiles = []
        
        for idx, smiles in enumerate(df['smiles']):
            if pd.isna(smiles):
                invalid_indices.append(idx)
                invalid_smiles.append('(missing)')
                continue
            
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                invalid_indices.append(idx)
                invalid_smiles.append(str(smiles))
        
        total = len(df)
        invalid_count = len(invalid_indices)
        valid_count = total - invalid_count
        
        self.log_message(f"SMILES validation: {valid_count}/{total} valid, {invalid_count} invalid")
        
        if invalid_count > 0:
            # Show first few invalid SMILES as examples
            examples = invalid_smiles[:5]
            example_str = ', '.join([f"'{s}'" for s in examples])
            if invalid_count > 5:
                example_str += f" ... and {invalid_count - 5} more"
            
            self.log_message(f"Invalid SMILES examples: {example_str}")
            
            # Ask user what to do
            response = messagebox.askyesno(
                "Invalid SMILES Detected",
                f"Found {invalid_count} invalid SMILES out of {total} entries.\n\n"
                f"Examples: {example_str}\n\n"
                f"Do you want to continue anyway?\n"
                f"(Invalid entries will be removed from processing)"
            )
            
            if not response:
                return None
            
            # Drop rows with invalid SMILES
            df = df.drop(df.index[invalid_indices]).reset_index(drop=True)
            self.log_message(f"Dropped {invalid_count} rows with invalid SMILES")
            self.log_message(f"Remaining rows: {len(df)}")
        
        return df
    
    def create_widgets(self):
        # Main container with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights for resizing
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        # Title
        title_label = ttk.Label(main_frame, text="Spectral Denoising Frontend", 
                               font=('Helvetica', 16, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=10)
        
        # File Input Section
        input_frame = ttk.LabelFrame(main_frame, text="Input File", padding="10")
        input_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        input_frame.columnconfigure(1, weight=1)
        
        ttk.Label(input_frame, text="File Type:").grid(row=0, column=0, sticky=tk.W, pady=5)
        file_type_frame = ttk.Frame(input_frame)
        file_type_frame.grid(row=0, column=1, sticky=tk.W, pady=5)
        ttk.Radiobutton(file_type_frame, text="MSP", variable=self.file_type, 
                       value="msp").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(file_type_frame, text="CSV", variable=self.file_type, 
                       value="csv").pack(side=tk.LEFT, padx=5)
        
        ttk.Label(input_frame, text="File Path:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(input_frame, textvariable=self.input_file_path, width=50).grid(
            row=1, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        ttk.Button(input_frame, text="Browse...", command=self.browse_input_file).grid(
            row=1, column=2, pady=5)
        
        ttk.Button(input_frame, text="Load File", command=self.load_file).grid(
            row=2, column=1, pady=10)
        
        # Parameters Section
        param_frame = ttk.LabelFrame(main_frame, text="Processing Parameters", padding="10")
        param_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(param_frame, text="Mass Tolerance (Da):").grid(row=0, column=0, sticky=tk.W, pady=5)

        tol_frame = ttk.Frame(param_frame)
        tol_frame.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)

        ttk.Radiobutton(tol_frame, text="Orbitrap (0.005)", variable=self.mass_mode,
                       value="orbitrap", command=self.on_mass_mode_change).grid(row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(tol_frame, text="ToF (0.01)", variable=self.mass_mode,
                       value="tof", command=self.on_mass_mode_change).grid(row=0, column=1, sticky=tk.W, padx=(10, 0))
        ttk.Radiobutton(tol_frame, text="Custom:", variable=self.mass_mode,
                       value="custom", command=self.on_mass_mode_change).grid(row=0, column=2, sticky=tk.W, padx=(10, 0))

        self.custom_tol_entry = ttk.Entry(tol_frame, textvariable=self.mass_tolerance, width=10, state='disabled')
        self.custom_tol_entry.grid(row=0, column=3, sticky=tk.W, padx=(5, 0))
        
        # Add more parameters here as needed
        # Example: ttk.Label(param_frame, text="Other Param:").grid(row=1, column=0, sticky=tk.W, pady=5)
        
        # Output Section
        output_frame = ttk.LabelFrame(main_frame, text="Output File", padding="10")
        output_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        output_frame.columnconfigure(1, weight=1)
        
        ttk.Label(output_frame, text="Save Path:").grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Entry(output_frame, textvariable=self.output_file_path, width=50).grid(
            row=0, column=1, sticky=(tk.W, tk.E), pady=5, padx=5)
        ttk.Button(output_frame, text="Browse...", command=self.browse_output_file).grid(
            row=0, column=2, pady=5)
        
        # Action Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=15)
        
        ttk.Button(button_frame, text="Process Data", command=self.process_data, 
                  style='Accent.TButton').pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Exit", command=self.root.quit).pack(side=tk.LEFT, padx=5)
        
        # Progress Section
        progress_frame = ttk.LabelFrame(main_frame, text="Progress", padding="10")
        progress_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.progress_bar = ttk.Progressbar(progress_frame, mode='determinate', maximum=100)
        self.progress_bar.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5)
        
        self.progress_label = ttk.Label(progress_frame, text="0%")
        self.progress_label.grid(row=0, column=1, padx=5)
        
        progress_frame.columnconfigure(0, weight=1)
        
        # Log/Output Section
        log_frame = ttk.LabelFrame(main_frame, text="Log", padding="10")
        log_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        main_frame.rowconfigure(6, weight=1)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        
        self.log_text = scrolledtext.ScrolledText(log_frame, height=10, width=80, 
                                                   state='disabled', wrap=tk.WORD)
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Initialize tolerance control state
        self.on_mass_mode_change()
    
    def log_message(self, message):
        """Add message to log window"""
        self.log_text.config(state='normal')
        self.log_text.insert(tk.END, message + '\n')
        self.log_text.see(tk.END)
        self.log_text.config(state='disabled')
        self.root.update_idletasks()
    
    def update_progress(self, current, total):
        """Update progress bar with current progress"""
        if total > 0:
            percentage = (current / total) * 100
            self.progress_bar['value'] = percentage
            self.progress_label.config(text=f"{int(percentage)}% ({current}/{total})")
            self.root.update_idletasks()

    def on_mass_mode_change(self):
        """Handle mass tolerance preset selection"""
        mode = self.mass_mode.get()
        if mode == 'orbitrap':
            self.mass_tolerance.set(0.005)
            self.custom_tol_entry.config(state='disabled')
            self.log_message(f"Mass tolerance set to: 0.005 Da (Orbitrap)")
        elif mode == 'tof':
            self.mass_tolerance.set(0.01)
            self.custom_tol_entry.config(state='disabled')
            self.log_message(f"Mass tolerance set to: 0.01 Da (ToF)")
        else:
            # Custom entry enabled
            self.custom_tol_entry.config(state='normal')
            self.custom_tol_entry.focus_set()
            self.log_message(f"Mass tolerance mode: Custom (current value: {self.mass_tolerance.get()} Da)")
    
    def browse_input_file(self):
        """Open file dialog to select input file"""
        file_types = [
            ("MSP files", "*.msp"),
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
        filename = filedialog.askopenfilename(
            title="Select Input File",
            filetypes=file_types
        )
        if filename:
            self.input_file_path.set(filename)
            # Auto-detect file type
            if filename.endswith('.msp'):
                self.file_type.set('msp')
            elif filename.endswith('.csv'):
                self.file_type.set('csv')
    
    def browse_output_file(self):
        """Open file dialog to select output file location"""
        file_types = [
            ("MSP files", "*.msp"),
            ("CSV files", "*.csv"),
            ("All files", "*.*")
        ]
        filename = filedialog.asksaveasfilename(
            title="Save Output File",
            filetypes=file_types,
            defaultextension=".csv"
        )
        if filename:
            self.output_file_path.set(filename)
    
    def load_file(self):
        """Load the input file"""
        file_path = self.input_file_path.get()
        
        if not file_path:
            messagebox.showerror("Error", "Please select an input file")
            return
        
        if not os.path.exists(file_path):
            messagebox.showerror("Error", "File does not exist")
            return
        
        try:
            self.log_message(f"Loading file: {file_path}")
            
            if self.file_type.get() == 'msp':
                self.loaded_data = file_io.read_msp(file_path)
                self.log_message(f"Successfully loaded MSP file with {len(self.loaded_data)} spectra")
            elif self.file_type.get() == 'csv':
                self.loaded_data = file_io.read_df(file_path)
                self.log_message(f"Successfully loaded CSV file with {len(self.loaded_data)} rows")
            
            # Standardize column names
            if isinstance(self.loaded_data, pd.DataFrame):
                self.loaded_data = self.standardize_columns(self.loaded_data)
                
                # Try to extract SMILES from comments if missing
                self.loaded_data = self.fill_missing_smiles(self.loaded_data)
                
                # Validate SMILES structures
                self.loaded_data = self.validate_smiles(self.loaded_data)
                if self.loaded_data is None:
                    self.log_message("Data loading cancelled by user due to invalid SMILES")
                    return
                
                # Validate required columns
                required_columns = ['peaks', 'smiles', 'adducts']
                missing_columns = [col for col in required_columns if col not in self.loaded_data.columns]
                
                if missing_columns:
                    error_msg = f"Missing required columns: {', '.join(missing_columns)}\n\n"
                    error_msg += f"Available columns: {', '.join(self.loaded_data.columns.tolist())}\n\n"
                    error_msg += "Please select a file with 'peaks', 'smiles', and 'adducts' columns."
                    self.log_message(f"Validation failed: Missing columns {missing_columns}")
                    self.loaded_data = None
                    messagebox.showerror("Invalid File Format", error_msg)
                    return
                
                # Show total number of spectra loaded
                self.log_message(f"Total spectra loaded: {len(self.loaded_data)}")
            
            messagebox.showinfo("Success", f"File loaded!\n Total spectra loaded: {len(self.loaded_data)}")
            
        except Exception as e:
            self.log_message(f"Error loading file: {str(e)}")
            messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
    
            
    def process_data(self):
        """Process the loaded data"""
        if self.loaded_data is None:
            messagebox.showerror("Error", "Please load a file first")
            return
        
        if not self.output_file_path.get():
            messagebox.showwarning("Warning", "No output path specified")
            return
        
        # Run processing in a separate thread to avoid freezing GUI
        thread = threading.Thread(target=self._process_data_thread)
        thread.start()
    
    def _process_data_thread(self):
        """Background thread for processing data"""
        try:
            self.progress_bar['value'] = 0
            self.progress_label.config(text="0%")
            self.log_message("Starting data processing...")
            
            # Example processing - modify this based on your actual needs
            # This is where you'd call your spectral_denoising functions
            
            mass_tol = self.mass_tolerance.get()
            self.log_message(f"Using mass tolerance: {mass_tol}")
            
            # Example: Process the data with progress callback
            # If your data has msms, smiles, and adduct columns:
            # self.results = denoising_pipeline(self.loaded_data, mass_tol, progress_callback=self.update_progress)
            self.results = denoising_pipeline(self.loaded_data, mass_tol)
            # For now, just demonstrate saving
            output_path = self.output_file_path.get()
            
            if isinstance(self.loaded_data, pd.DataFrame):
                if output_path.endswith('.csv'):
                    # self.loaded_data.to_csv(output_path, index=False)
                    file_io.save_df(self.results, output_path)
                    self.log_message(f"Saved results to: {output_path}")
                else:
                    file_io.write_to_msp(self.results, output_path, msms_col = 'denoised_peaks')# make sure to change this column!!!
            
            self.progress_bar['value'] = 100
            self.progress_label.config(text="100%")
            self.log_message("Processing completed successfully!")
            messagebox.showinfo("Success", "Data processing completed!")
            
        except Exception as e:
            self.progress_bar['value'] = 0
            self.progress_label.config(text="Error")
            self.log_message(f"Error during processing: {str(e)}")
            messagebox.showerror("Error", f"Processing failed:\n{str(e)}")
    
    def clear_all(self):
        """Clear all inputs and outputs"""
        self.input_file_path.set("")
        self.output_file_path.set("")
        self.loaded_data = None
        self.log_text.config(state='normal')
        self.log_text.delete(1.0, tk.END)
        self.log_text.config(state='disabled')
        self.log_message("Cleared all data")


def main():
    root = tk.Tk()
    app = SpectralDenoisingGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
