import numpy as np
import sys, os
sys.path.append(r'C:\Users\User\Documents\03_software\python')
from PyMS import Sequences as Seq
from PyMS import CalculateMW as calcMW


class spectrum(object):
    def __init__(self, start_nm=190, stop_nm=260, step_nm=1):
        """
        Initialise the class

        Args:
            start_nm: Value to start the spectrum
            stop_nm: Value to end the spectrum
            step_nm: Step-width in nm
        """

        # init the wavelength dimension
        self.wavelength_nm = np.arange(start_nm, stop_nm + 1, step_nm)
        # init the ellipticity as nan
        self.ellipticity_mdeg = np.full(1+int((stop_nm - start_nm)/step_nm), np.nan)
        # init the detector coltage as nan
        self.detector_V = np.full(1+int((stop_nm - start_nm)/step_nm), np.nan)
        # init mean-residue weighted ellipticity
        self.mrwe_deg_cm_2_dmol = None

    def load_jasco(self, filepath):
        """
        Open a JASCO results file and store the values in self

        Args:
            infile: path to infile
        Returns:
            self.wavelength_nm (np.array)
            self.CD_mdeg (np.array)
            self.HT_V (np.array)
        """

        with open(filepath,'r') as inf:
            line = inf.readline().strip()
            datablock = False
            while line != '':
                if line.startswith('XYDATA'):
                    datablock = True
                elif datablock is True:
                    x, y1, y2 = [i.strip() for i in line.strip().split()]
                    for idx, wl_self in enumerate(self.wavelength_nm):
                        if float(x) == float(wl_self):
                            self.ellipticity_mdeg[idx] = y1
                            self.detector_V[idx] = y2
                            break

                line = inf.readline().strip()

    def filter_by_detector(self, minV=0, maxV=1000, inplace=True):
        """
        Remove spectral data by including only those measurements lying within
        a range of detector voltages
        
        Args:
            minV: minimum detector voltage (default 0)
            maxV: maximum detector volatge (default 1000)
        """
        
        nan_idx  = (self.detector_V > maxV) | (self.detector_V < minV)
        
        a = self.ellipticity_mdeg
        a[nan_idx] = np.nan
        b = None
        
        if self.mrwe_deg_cm_2_dmol is not None:
            b = self.mrwe_deg_cm_2_dmol
            b[nan_idx] = np.nan
        
        if inplace is True:
            self.ellipticity_mdeg, self.mrwe_deg_cm_2_dmol = a, b
        else:
            return a,b

    def get_MRWE(self, MRW_g_per_mol, pathlength_cm, conc_mg_per_mL):
        """
        Calculate mean residue-weighted ellipticity from the raw-data

        Args:
            MRW_g_per_mol: Mass of the protein divided by the amino acid count N - 1
            pathlength_cm: cuvette path length in cm
            conc_mg_per_mL: protein concentration in mg/mL
        Returns:
            self.MRWE
        """
        
        print(MRW_g_per_mol)
        print(conc_mg_per_mL)
        print(pathlength_cm)
        # Calculate mean residue weighted ellipticity
        self.mrwe_deg_cm_2_dmol = self.ellipticity_mdeg * (142.5*3300/5610051) * MRW_g_per_mol / (conc_mg_per_mL * pathlength_cm)

    def add(self, new_spectrum):
        if np.array_equal(self.wavelength_nm, new_spectrum.wavelength_nm):
            self.ellipticity_mdeg = self.ellipticity_mdeg + new_spectrum.ellipticity_mdeg
            if self.mrwe_deg_cm_2_dmol is not None and new_spectrum.mrwe_deg_cm_2_dmol is not None:
                self.mrwe_deg_cm_2_dmol = self.mrwe_deg_cm_2_dmol + new_spectrum.mrwe_deg_cm_2_dmol
        else:
            raise Exception('the wavelength ranges and spacings of the two spectra must match')

    def substract (self, new_spectrum, inplace=True):
        if np.array_equal(self.wavelength_nm, new_spectrum.wavelength_nm):
            if inplace is True:
                self.ellipticity_mdeg = self.ellipticity_mdeg - new_spectrum.ellipticity_mdeg
                if self.mrwe_deg_cm_2_dmol is not None and new_spectrum.mrwe_deg_cm_2_dmol is not None:
                    self.mrwe_deg_cm_2_dmol = self.mrwe_deg_cm_2_dmol - new_spectrum.mrwe_deg_cm_2_dmol
            else:
                a = self.ellipticity_mdeg - new_spectrum.ellipticity_mdeg
                b = None
                if self.mrwe_deg_cm_2_dmol is not None and new_spectrum.mrwe_deg_cm_2_dmol is not None:
                    b = self.mrwe_deg_cm_2_dmol - new_spectrum.mrwe_deg_cm_2_dmol
                return a, b

        else:
            raise Exception('the wavelength ranges and spacings of the two spectra must match')

# Usage example
#myspectrum = spectrum()
#mybaseline = spectrum()
#
#myspectrum.load_jasco(r'C:\Users\User\Downloads\cd\20191029_jb_syt1_c2astar_go_buffer.txt')
#mybaseline.load_jasco(r'C:\Users\User\Downloads\cd\20191030_jb_baseline_buffer.txt')
#
## Get the protein sequence from the fasta-file
#IDList, sequenceList = Seq.ReadFasta('syt1_soluble.fasta')
#for idx, identifier in enumerate(IDList):
#    if identifier == 'Syt1 C2A GO':
#        protein_sequence = sequenceList[idx]
#
## Calculate the mean residue mass of the protein from the sequence
#mean_residue_mass = calcMW.mass_from_sequence(protein_sequence, mono=False, returnAverage=False) / (len(protein_sequence) -1)
#
#myspectrum.get_MRWE(mean_residue_mass, 0.01, 0.935)
#myspectrum.substract(mybaseline)