import MDAnalysis as mda
from MDAnalysis import analysis
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.contacts import contact_matrix
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
from multiprocessing import Pool
from multiprocessing import cpu_count
import random
import seaborn as sns
import nglview as nv
# from numba import jit, cuda
# from timeit import default_timer as timer  

quotes = ['"Lets-a go" - Chris Pratt','"I look very much forward to showing my financials, because they are huge." - Donald Trump','"I did not copy Gromacs" - Some rad guy, probably']

def help(self,func=''):
    if func == '':
        print('Specify a function for help.')
        print('e.g. heatmapper.help("frame_to_time")')
    if func == 'frame_to_time':
        print('---Arguments---')
        print('frame_in - The simulation frame for which you want to know the coresponding time')
    if func == 'time_to_frame':
        print('---Arguments---')
        print('time_in - The simulation time for which you want to know the coresponding frame')
    if func == 'cont_per_frame':
        print('---Arguments---')
        print('- This function is mainly for the multi-threading function.')
        print('frame_index - Frame to analyse')
        print('cont_dist - Cutoff distance to define a contact')
        print('segid - Segid to analyse. Contacts are far easier to keep track of if you treat each segid sequentially.')
        print('probe_sel - Selection for the probe moleule to count contact from. This uses MDA selection language.')
        print('carbon - Boolean. If False then carbons are not included for contacts. Can speed up very large systems.')
    if func == 'cont_pro':
        print('---Arguments---')
        print('cont_dist - Cutoff distance to define a contact. cont_dist=3.3')
        print('probe_sel - Selection for the probe moleule to count contact from. This uses MDA selection language. probe_sel="protein"')
        print('carbon - Boolean. If False then carbons are not included for contacts. Can speed up very large systems.carbon=True')
        print('start - Frame to start contact analysis at. start=0')
        print('stop - Frame to stop contact analysis at. stop=-1')
        print('skip - Frames to skip. Can speed up for large systems or long trajectories. skip=1')

class traj_analysis:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,pdb_filename,xtc_filename):
        self.pdb_filename = pdb_filename
        self.xtc_filename = xtc_filename
        self.mda_universe = mda.Universe(pdb_filename, xtc_filename) # <--- Load in universe
        # self.mda_universe_coords = mda.Universe(pdb_filename)
        self.protein = self.mda_universe.select_atoms('protein') # <--- Find the protein
        self.n_jobs = cpu_count() # <--- CPUs available for multithreading
        self.n_frames_tot = self.mda_universe.trajectory.trajectory.n_frames # <--- Find total number of frames in the trajectory
        rnd_qt = random.randint(0,int(len(quotes)-1))
        print("- Run the heatmapper.help() function if you're lost")
        print(str(quotes[rnd_qt]))
    def frame_to_time(self,frame_in):
        self.time = frame_in*self.mda_universe.trajectory.dt
        return self.time
    def time_to_frame(self,time_in):
        self.frame = time_in/self.mda_universe.trajectory.dt
        return self.frame
    def cont_per_frame(self,frame_index,cont_dist,segid,analyte,probe_sel,carbon): # <--- The actual function which is executed on each CPU
        print('Reached frame '+str(frame_index)+'. Final frame is '+str(self.frame_values[-1])+'.           ',end='\r')
        self.mda_universe.trajectory[frame_index] # <--- Select the frame for analysis
        residue_contacts=[] # <--- Create an empty array to store the contacts per residue in
        for resid_iter in self.pro_resid: # <--- Iterate throught resids
            if carbon == True:
                group_A = self.mda_universe.select_atoms('('+analyte+') and (resid '+str(resid_iter)+')') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(5*cont_dist)+' resid '+str(resid_iter)+') and '+probe_sel) # <--- Select the analyte around the resid
            elif carbon == False:
                group_A = self.mda_universe.select_atoms('('+analyte+') and (resid '+str(resid_iter)+') and (not name C)') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(5*cont_dist)+' resid '+str(resid_iter)+') and (not name C) and '+probe_sel) # <--- Select the analyte around the resid
            distances = mda.analysis.distances.distance_array(group_A.positions, group_B.positions) # <--- Get distances
            contact_count = np.count_nonzero(distances <= cont_dist) # <--- Count the number of distances under the cutoff
            residue_contacts.append(contact_count) # <--- Add the number of contacts for that residus
        self.cont_frames_done.append(frame_index)
        return residue_contacts
    def cont_pro(self,cont_dist=3.3,analyte='protein',probe_sel='protein',carbon=True,start=0,stop=-1,skip=1):
        self.cont_frames_done = []
        self.raw_cont_data = [] # <--- Create empty array for the chain data
        segid_iter = 0
        segid_iter = 0
        all_dfs = []
        segids_in_pro = self.mda_universe.select_atoms(analyte).segids
        print('Selected analyte')
        self.pro_segid = []
        for segid in segids_in_pro:
            if segid not in self.pro_segid:
                self.pro_segid.append(segid) # <--- Get a list of segids in the protein
        print('Got segids')
        resids_in_pro = self.mda_universe.select_atoms(analyte).resids
        self.pro_resid = []
        for resid in resids_in_pro:
            if resid not in self.pro_resid:
                self.pro_resid.append(resid) # <--- Get a list of resids in the protein
        print('Got resids')
        for segid in self.pro_segid: # <--- Iterate through chains
            print('Chain: '+str(segid),end='\n') # <--- Update console
            run_per_frame = partial(self.cont_per_frame, cont_dist=cont_dist,segid=segid,analyte=analyte,probe_sel=probe_sel,carbon=carbon) # <--- Set the per frame function
            frame_values = np.arange(self.mda_universe.trajectory.n_frames) # <--- Get all frames
            frame_values = frame_values[start:stop:skip] # <--- Select every nth frame for reducing run time
            self.frame_values = frame_values
            print(str(len(frame_values))+' frames to analyse on '+str(self.n_jobs)+' cores.',end='\n')
            print('First frame: '+str(self.frame_values[0])+' - Final frame: '+str(self.frame_values[-1]))
            with Pool(self.n_jobs) as worker_pool: # <--- Create a pool of CPUs to use
                print('Running analysis')
                result = worker_pool.map(run_per_frame, frame_values) # <--- Run the per frame function on a select CPU
                print('Analysis done')
            self.raw_cont_data.append(np.asarray(result).T) # <--- Transpose the results into a single array
            df = pd.DataFrame({'resids':self.pro_resid})
            time_data_dfs = []
            for iteration in range(len(self.raw_cont_data[0][0])):
                time_data = []
                for resid_iter in range(len(self.raw_cont_data[0])):
                    time_data.append(self.raw_cont_data[segid_iter][resid_iter][iteration])
                time_data_dfs.append(pd.DataFrame({'frame_'+str(frame_values[iteration]):time_data}))
            df_concat = pd.concat([df, *time_data_dfs], axis=1)
            all_dfs.append(df_concat)
            segid_iter+=1
        return all_dfs
    def gen_map(self,df_data,res_op=0.6,res_rep='ball+stick'):
        segid_iter = 0
        self.heatmaps = []
        for segid in self.pro_segid:
            df = df_data[segid_iter].copy()
            conts_per_resid = []
            for i in range(len(df.resids)):
                conts_per_resid.append(df.iloc[i,1:].sum())
            self.min_conts = min(conts_per_resid)
            self.max_conts = max(conts_per_resid)
            conts_per_resid_minus_min = conts_per_resid-self.min_conts # <--- Scale starts at 0
            conts_per_resid_std = ((conts_per_resid_minus_min/max(conts_per_resid_minus_min))*100).astype(int) # <--- 100 colour values
            sequential_colors = sns.color_palette("flare", 101).as_hex()
            segid_sel = 'protein and segid '+str(segid)
            segid_act = self.mda_universe.select_atoms(segid_sel)
            view = nv.show_mdanalysis(segid_act,default_representation=False)
            view.add_representation('cartoon', selection='all', color='grey', opacity=1)
            resid_iter = 0
            for resid in df.resids:
                view.add_representation(res_rep, selection=str(resid), color=sequential_colors[conts_per_resid_std[resid_iter]], opacity=res_op)
                resid_iter+=1
            view._view_height = '600px'
            self.heatmaps.append(view)
        return self.heatmaps, sequential_colors
    def res_cont_freq(self,df_data,segid_identifier='protein'):
        AA_contacts_dict = {
            'CYS': 0, 'ASP': 0, 'SER': 0, 'GLN': 0, 
            'LYS': 0, 'ILE': 0, 'PRO': 0, 'THR': 0, 
            'PHE': 0, 'ASN': 0, 'GLY': 0, 'HIS': 0, 
            'LEU': 0, 'ARG': 0, 'TRP': 0, 'ALA': 0, 
            'VAL':0, 'GLU': 0, 'TYR': 0, 'MET': 0}
        for resid in self.pro_resid:
            resname = (self.mda_universe.select_atoms('resid '+str(resid)).residues.resnames[0]).upper()
            AA_contacts_dict[resname] += df_data.query('resids == @resid').iloc[0,1:].values.sum()
        dic_out = {}
        for x,y in AA_contacts_dict.items():
            if y != 0:
                dic_out[x] = y # <--- Remove entries with 0 contacts. Very likely a residue which is not in the protein.
        for resname in dic_out.keys():
            occourance_count = 0
            for item in self.mda_universe.select_atoms(segid_identifier).residues.resnames:
                if item.upper() == resname.upper():
                    occourance_count += 1
            if occourance_count == 0:
                dic_out[resname] = 0
            else:
                dic_out[resname] = dic_out[resname] / occourance_count
        self.AA_conts = pd.DataFrame.from_dict(dic_out, orient='index', columns=['contacts'])
        return self.AA_conts