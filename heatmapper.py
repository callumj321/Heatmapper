import MDAnalysis as mda
from MDAnalysis import analysis
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.contacts import contact_matrix
import mdtraj as md
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
from multiprocessing import Pool
from multiprocessing import cpu_count
import random
import seaborn as sns
import nglview as nv
import time
# from numba import jit, cuda
# from timeit import default_timer as timer  

quotes = ['"Lets-a go" - Chris Pratt','"I look very much forward to showing my financials, because they are huge." - Donald Trump','"I did not copy Gromacs" - Some rad guy, probably',
'"ff at 15" - Abraham Lincoln']

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
        print('---Usage---')
        print('This function will take an input coordinate and trajectory file and put out a summary of the contacts from a probe molecule. By default the analyte and probe selections are "protein", so will return inter-protein contacts.')
        print('---Arguments---')
        print('cont_dist - Cutoff distance to define a contact. cont_dist=3.3')
        print('probe_sel - Selection for the probe moleule to count contact from. This uses MDA selection language. probe_sel="protein"')
        print('carbon - Boolean. If False then carbons are not included for contacts. Can speed up very large systems.carbon=True')
        print('start - Frame to start contact analysis at. start=0')
        print('stop - Frame to stop contact analysis at. stop=-1')
        print('skip - Frames to skip. Can speed up for large systems or long trajectories. skip=1')
        print('---Output---')
        print('all_dfs - A list containing a contact dataframe for each protein chain in the input structure.')

class traj_analysis:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,pdb_filename,xtc_filename,analyte='protein'):
        self.pdb_filename = pdb_filename
        self.xtc_filename = xtc_filename
        print('Loading in mda universe...')
        self.mda_universe = mda.Universe(pdb_filename,xtc_filename) # <--- Load in universe
        print('Done.')
        # self.mdtraj_universe = md.load(pdb_filename,xtc_filename)
        # self.mda_universe_coords = mda.Universe(pdb_filename)
        self.analyte_loaded = self.mda_universe.select_atoms(analyte)
        self.analyte_resids = np.unique(self.analyte_loaded.residues.resids)
        self.n_jobs = cpu_count() # <--- CPUs available for multithreading
        self.n_frames_tot = self.mda_universe.trajectory.trajectory.n_frames # <--- Find total number of frames in the trajectory
        rnd_qt = random.randint(0,int(len(quotes)-1))
        print(str(quotes[rnd_qt])+'\n')
        print(str(len(np.unique(self.analyte_loaded.segids)))+' segids loaded.')
        print(str(len(self.analyte_resids))+' resids per segid loaded.\n')
        print("\n- Run the heatmapper.help() function if you're lost")
    def frame_to_time(self,frame_in):
        self.time = frame_in*self.mda_universe.trajectory.dt
        return self.time
    def time_to_frame(self,time_in):
        self.frame = time_in/self.mda_universe.trajectory.dt
        return self.frame
    def pre_cont_pro(self,start,stop,skip):
        frame_values = np.arange(self.mda_universe.trajectory.n_frames) # <--- Get all frames
        print('Frames: '+str(int(len(frame_values[start:stop:skip]))))
        print('Cores: '+str(int(self.n_jobs)))
        print('Passes: '+str(int(len(frame_values[start:stop:skip])/self.n_jobs)))
    def cont_per_frame(self,frame_index,cont_dist,segid,analyte,probe_sel,carbon): # <--- The actual function which is executed on each CPU
        start_time = time.time()
        print('Reached frame '+str(frame_index)+'           ',end='\r')
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
        self.cont_pro_run_times.append(time.time() - start_time)
        return residue_contacts
    def cont_pro(self,cont_dist=3.3,analyte='protein',probe_sel='protein',carbon=True,start=0,stop=-1,skip=1):
        self.cont_pro_run_times = []
        self.cont_frames_done = []
        self.raw_cont_data = [] # <--- Create empty array for the chain data
        segid_iter = 0
        segid_iter = 0
        self.contact_dfs = []
        segids_in_pro = self.mda_universe.select_atoms(analyte).segids
        print('Selected analyte')
        self.pro_segid = []
        for segid in segids_in_pro:
            if segid not in self.pro_segid:
                self.pro_segid.append(segid) # <--- Get a list of segids in the protein
        print('Got segids. ',str(len(self.pro_segid)),' segids to analyse.')
        resids_in_pro = self.mda_universe.select_atoms(analyte).resids
        self.pro_resid = []
        for resid in self.analyte_resids:
            if resid not in self.pro_resid:
                self.pro_resid.append(resid) # <--- Get a list of resids in the protein
        print('Got resids. ',str(len(self.pro_resid)),' resids to analyse.')
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
            self.contact_dfs.append(df_concat)
            segid_iter+=1
        # self.gen_map(self.contact_dfs,analyte)
        return self.contact_dfs
    def gen_map(self,df_data,analyte,res_op=1,res_rep='spacefill',res_rad=0.5,BB_only=True):
        # Standardise df so that the resids are the index.
        print('Standardising input')
        if 'resids' in df_data.columns.values:
            df_data = df_data.set_index('resids')
        print('Getting contact totals')
        conts_per_resid = []
        for resid in self.analyte_resids:
            conts_per_resid.append(df_data.query('index == @resid').values.sum())
        # return conts_per_resid
        print('Got contact totals',end='\r')
        self.min_conts = min(conts_per_resid)
        self.max_conts = max(conts_per_resid)
        conts_per_resid_minus_min = conts_per_resid-self.min_conts # <--- Scale starts at 0
        conts_per_resid_std = ((conts_per_resid_minus_min/max(conts_per_resid_minus_min))*100).astype(int) # <--- 100 colour values
        sequential_colors = sns.color_palette("flare", 101).as_hex()
        act = self.mda_universe.select_atoms(analyte)
        self.heatmap = nv.show_mdanalysis(act,default_representation=False)
        self.heatmap.add_representation('cartoon', selection='all', color='grey', opacity=1)
        resid_iter = 0
        if BB_only == True:
            for resid in self.analyte_resids:
                self.heatmap.add_representation(res_rep,selection='backbone and '+str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                resid_iter+=1
        else:
            for resid in self.analyte_resids:
                self.heatmap.add_representation(res_rep,selection=str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                resid_iter+=1
        self.heatmap._view_height = '800px'
        print('Done. Run \nobj.heatmap \nin jupyter to view. Run \nnv.write_html("index.html",[obj.heatmap]) \nTo save as an html you can open later.')
        return self.heatmap
    def AA_cont_freq(self,df_data,segid_identifier='protein'):
        print('Running')
        AA_contacts_dict = {
            'CYS': 0, 'ASP': 0, 'SER': 0, 'GLN': 0, 
            'LYS': 0, 'ILE': 0, 'PRO': 0, 'THR': 0, 
            'PHE': 0, 'ASN': 0, 'GLY': 0, 'HSD': 0, 
            'LEU': 0, 'ARG': 0, 'TRP': 0, 'ALA': 0, 
            'VAL':0, 'GLU': 0, 'TYR': 0, 'MET': 0}
        AA_het_at_dict = AA_contacts_dict.copy()
        AA_CH_ratio_dict = AA_contacts_dict.copy()
        AA_prop_dict = {
            'CYS': 'Special', 'ASP': '-ve', 'SER': 'Polar', 'GLN': 'Polar', 
            'LYS': '+ve', 'ILE': 'Hydrophobic', 'PRO': 'Special', 'THR': 'Polar', 
            'PHE': 'Hydrophobic', 'ASN': 'Polar', 'GLY': 'Special', 'HSD': '+ve', 
            'LEU': 'Hydrophobic', 'ARG': '+ve', 'TRP': 'Hydrophobic', 'ALA': 'Hydrophobic', 
            'VAL': 'Hydrophobic', 'GLU': '-ve', 'TYR': 'Hydrophobic', 'MET': 'Hydrophobic'}
        # Standardise df so that the resids are the index.
        print('Standardising input')
        if 'resids' in df_data.columns.values:
            df_data = df_data.set_index('resids')
        print('Getting hetro atom counts')
        for AA in AA_het_at_dict.keys():
            at_list = self.mda_universe.select_atoms('protein and resname '+str(AA)).residues[0].atoms.types
            for at in at_list:
                if at != 'H':
                    AA_het_at_dict[AA]+=1
                AA_CH_ratio_dict[AA] = ((at_list=='C').sum()+(at_list=='H').sum())/len(at_list)
        print('Getting contacts')
        for resid in self.analyte_resids: # <--- VERY SLOW
            resname = (self.mda_universe.select_atoms('resid '+str(resid)).residues.resnames[0]).upper()
            AA_contacts_dict[resname] += df_data.query('index == @resid').values.sum()
        # Divide the numnber of contacts by the numnber of that residue in the protein so that they're not over represented.
        print('Standardising contacts by number of residue type')
        for resname in AA_contacts_dict.keys():
            occourance_count = 0
            for item in self.mda_universe.select_atoms(segid_identifier).residues.resnames:
                if item.upper() == resname.upper():
                    occourance_count += 1
            if occourance_count == 0:
                AA_contacts_dict[resname] = 0
            else:
                AA_contacts_dict[resname] = AA_contacts_dict[resname] / occourance_count
        # Report the number of hetroatoms in each amino acid so that size can be considered.
        print('Reporting')
        self.AA_conts = pd.DataFrame.from_dict(AA_contacts_dict, orient='index', columns=['contacts'])
        self.AA_conts['het_count'] = AA_het_at_dict.values()
        self.AA_conts['prop'] = AA_prop_dict.values()
        self.AA_conts['CH_ratio'] = AA_CH_ratio_dict.values()
        print('BG. Reported.')
        return self.AA_conts
    def res_cont_freq(self,df_data):
        print('Standardising input')
        if 'resids' in df_data.columns.values:
            df_data = df_data.set_index('resids')
        resid_tot_conts = np.zeros(len(df_data.index.values))
        resid_iter = 0
        for resid in self.analyte_resids:
            resid_tot_conts[resid_iter] = df_data.query('index == @resid').values.sum()
            resid_iter += 1
        self.cont_per_res = pd.DataFrame({'resids':df_data.resids.values,'contacts':resid_tot_conts})
        return self.cont_per_res