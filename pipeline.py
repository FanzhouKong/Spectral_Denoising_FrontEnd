from spectral_denoising.spectral_denoising import spectral_denoising
from tqdm import tqdm
import multiprocessing as mp
from spectral_denoising.spectral_denoising import spectral_denoising_batch



def denoising_pipeline(input_df, mass_error):
    input_df.columns = [col.lower() for col in input_df.columns]
    peaks = input_df['peaks']
    smiles = input_df['smiles']
    adducts = input_df['adducts']
        
    denoised_peaks = spectral_denoising_batch(peaks, smiles, adducts, mass_tolerance=mass_error)
    
    output_df = input_df.copy()
    output_df['denoised_peaks'] = denoised_peaks
    
    return output_df