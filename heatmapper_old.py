import MDAnalysis as mda
from MDAnalysis import analysis
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.contacts import contact_matrix
# import mdtraj as md
import pandas as pd
import numpy as np
import multiprocessing
from multiprocessing import Pool
from multiprocessing import cpu_count
from functools import partial
import random
import seaborn as sns
import nglview as nv
import time
from datetime import datetime
from datetime import timedelta
# from numba import jit, cuda
# from timeit import default_timer as timer  

def help(func=''):
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
    def __init__(self,pdb_filename,xtc_filename,analyte='protein',probe='protein'):
        self.pdb_filename = pdb_filename
        self.xtc_filename = xtc_filename
        print('Loading in mda universe...')
        self.mda_universe = mda.Universe(pdb_filename,xtc_filename) # <--- Load in universe
        print('Done.')
        # self.mdtraj_universe = md.load(pdb_filename,xtc_filename)
        # self.mda_universe_coords = mda.Universe(pdb_filename)
        self.analyte_sel = analyte
        self.probe_sel = probe
        self.analyte_loaded = self.mda_universe.select_atoms(analyte)
        self.analyte_segids = np.unique(self.analyte_loaded.residues.segids)
        self.analyte_resids = np.unique(self.analyte_loaded.residues.resids)
        self.n_jobs = cpu_count() # <--- CPUs available for multithreading
        self.n_frames_tot = self.mda_universe.trajectory.trajectory.n_frames # <--- Find total number of frames in the trajectory
        self.run_pass = 1
        rnd_qt = random.randint(0,int(len(quotes)-1))
        # print('\n'+str(quotes[rnd_qt])+'\n')
        print(str(len(np.unique(self.analyte_loaded.segids)))+' segids loaded.')
        print(str(len(self.analyte_resids))+' resids per segid loaded.')
        print("\n- Run the heatmapper.help() function if you're lost")
    def frame_to_time(self,frame_in):
        self.time = frame_in*self.mda_universe.trajectory.dt
        return self.time
    def time_to_frame(self,time_in):
        self.frame = time_in/self.mda_universe.trajectory.dt
        return self.frame
    def pre_cont_pro(self,start,stop,skip):
        print('Frames: '+str(int(len(self.mda_universe.trajectory[start:stop:skip]))))
        print('Cores: '+str(int(self.n_jobs)))
        print('Passes: '+str(int(len(self.mda_universe.trajectory[start:stop:skip])/self.n_jobs)))
        print('Chains: '+str(int(len(self.analyte_segids))))
        print('Residues per chain: '+str(int(self.analyte_resids)))
    def cont_per_frame(self,frame_index,cont_dist,carbon,segid='A'): # <--- The actual function which is executed on each CPU
        # print('Reached frame '+str(frame_index)+'           ',end='\r')
        print(str(frame_index),end='\r')
        self.mda_universe.trajectory[frame_index] # <--- Select the frame for analysis
        residue_contacts=[] # <--- Create an empty array to store the contacts per residue in
        for resid_iter in self.analyte_resids: # <--- Iterate throught resids
            if carbon == True:
                group_A = self.mda_universe.select_atoms('('+str(self.analyte_sel)+') and (resid '+str(resid_iter)+')') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(2*cont_dist)+' resid '+str(resid_iter)+') and ('+str(self.probe_sel)+')') # <--- Select the analyte around the resid
            elif carbon == False:
                group_A = self.mda_universe.select_atoms('('+str(self.analyte_sel)+') and (resid '+str(resid_iter)+') and (not name C)') # <--- Select the correct resid
                group_B = self.mda_universe.select_atoms('(around '+str(2*cont_dist)+' resid '+str(resid_iter)+') and ('+str(self.probe_sel)+') and (not name C)') # <--- Select the analyte around the resid
            distances = mda.analysis.distances.distance_array(group_A.positions, group_B.positions) # <--- Get distances
            contact_count = np.count_nonzero(distances <= cont_dist) # <--- Count the number of distances under the cutoff
            residue_contacts.append(contact_count) # <--- Add the number of contacts for that residus
        print('b')
        df_out = pd.read_csv('contact_analysis_data_'+str(segid)+'.csv')
        df_out[str(frame_index)] = residue_contacts
        print('Saving...')
        df_out.to_csv('contact_analysis_data_'+str(segid)+'.csv',index=False)
        print('Saved...)')
        self.run_pass += 1
    def cont_pro(self,cont_dist=3.3,carbon=True,start=0,stop=-1,skip=1):
        print(str(int(len(self.analyte_segids))),' segids to analyse.')
        print(str(int(len(self.analyte_resids))),' resids to analyse.')
        segid_iter = 0
        for segid in self.analyte_segids: # <--- Iterate through chains
            df_out = pd.DataFrame({'Resid':self.analyte_resids})
            df_out.to_csv('contact_analysis_data_'+str(segid)+'.csv',index=False)
            print('Chain: '+str(segid),end='\n') # <--- Update console
            run_per_frame = partial(self.cont_per_frame, cont_dist=cont_dist,carbon=carbon,segid=segid) # <--- Set the per frame function
            frame_values = np.arange(self.mda_universe.trajectory.n_frames) # <--- Get all frames
            self.analysis_frame_values = frame_values[start:stop:skip] # <--- Select every nth frame for reducing run time
            print(str(len(self.analysis_frame_values))+' frames to analyse on '+str(self.n_jobs)+' cores.',end='\n') # <--- Update console
            print('First frame: '+str(self.analysis_frame_values[0])+' - Final frame: '+str(self.analysis_frame_values[-1]))
            self.start_time = datetime.now().strftime("%H:%M:%S")
            print('Running analysis - Start time: ',str(self.start_time))
            print('Working on frame...')
            with Pool(self.n_jobs) as worker_pool: # <--- Create a pool of CPUs to use
                worker_pool.map(run_per_frame, self.analysis_frame_values) # <--- Run the per frame function on a select CPU
            stop_time = str(datetime.now().strftime("%H:%M:%S"))
            print('Analysis done - Stop time: ',stop_time)
            segid_iter+=1
            print('\n')
    def gen_map(self,df_data,res_op=1,res_rep='spacefill',res_rad=0,BB_only=True):
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
        act = self.mda_universe.select_atoms(self.analyte_sel)
        self.heatmap = nv.show_mdanalysis(act,default_representation=True)
        self.heatmap.add_representation('trace', selection='all', color='grey', opacity=1)
        resid_iter = 0
        if res_rad == 0:
            if BB_only == True:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection='backbone and '+str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op)
                    resid_iter+=1
            else:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection=str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op)
                    resid_iter+=1
        else:
            if BB_only == True:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection='backbone and '+str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                    resid_iter+=1
            else:
                for resid in self.analyte_resids:
                    self.heatmap.add_representation(res_rep,selection=str(resid),color=sequential_colors[conts_per_resid_std[resid_iter]],opacity=res_op,radius=res_rad)
                    resid_iter+=1
        self.heatmap._view_height = '400px'
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
    def display(self,selection='protein',spec_sel='all',res_op=1,res_rep='spacefill',res_rad=0,def_rep=True):
        act = self.mda_universe.select_atoms(str(selection))
        self.show_display = nv.show_mdanalysis(act,default_representation=def_rep)
        if res_rad == 0:
            self.show_display.add_representation(res_rep,selection=spec_sel,opacity=res_op)
        else:
            self.show_display.add_representation(res_rep,selection=spec_sel,opacity=res_op,radius=res_rad)
        self.show_display._view_height = '800px'
        self.show_display.center()
        return self.show_display