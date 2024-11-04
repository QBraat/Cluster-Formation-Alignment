
from cc3d.core.PySteppables import *
from pathlib import Path
import random
import os
import glob
from time import localtime

Np = int({{XDIM}}/10)
kappa = {{KAPPA}}
tau = {{TAU}}
frac = {{FRACTIONACTIVE}}
gMM = {{GMM}}
gEE = 0.0
gEM = 0.0




class ClusterCountingAlgorithmSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
     

    def start(self):
        global write_clusters, write_cells

        pixels = int(10)
        for i in range(Np):
            for j in range(Np):
                    self.cell_field[i*pixels:(i+1)*pixels, j*pixels:(j+1)*pixels,0] = self.new_cell(self.EPITHELIAL)
         
        Number_Mesenchymal = int(frac*Np*Np)
        
        while Number_Mesenchymal > 0: 
            cellid = random.randint(1,Np*Np) 
            for cell in self.cell_list:
                if cell.id == cellid and cell.type != self.MESENCHYMAL:
                    cell.type = self.MESENCHYMAL
                    Number_Mesenchymal -=1
                else: 
                    Number_Mesenchymal = Number_Mesenchymal
        
        volume_per_cell = int(self.dim.x*self.dim.y/(Np*Np))
        
        for cell in self.cell_list_by_type(self.EPITHELIAL,self.MESENCHYMAL):
            
            if cell.type == self.EPITHELIAL:
                cell.targetVolume = volume_per_cell
                cell.lambdaVolume = 2.0
            else:
                cell.targetVolume = volume_per_cell
                cell.lambdaVolume = 2.0
                
                
                
                
                
        write_clusters, file_path = self.open_file_in_simulation_output_folder('data_clusters_N%03d_kappa%04d_tau%04d_frac%03d_10000gMM%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_clusters is None:
            return   
           
        write_cells, file_path = self.open_file_in_simulation_output_folder('data_cells_N%03d_kappa%04d_tau%04d_frac%03d_10000gMM%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_cells is None:
            return    
             
            
        print('mcs cluster_id cluster_size',file=write_clusters)   
        print('mcs cell.id cell.type xCOM yCOM zCOM pol_angle cluster.id',file=write_cells)
                
                
                
                
        write_metadata, file_path = self.open_file_in_simulation_output_folder('metadata_N%03d_kappa%04d_tau%04d_frac%03d_10000gMM%05d.dat'%(int(Np),int(kappa),int(tau),int(frac*100),int(gMM*10000)), mode='w')
        if write_metadata is None:
            return   
            
        print('Manually set parameters  ',file=write_metadata)
        print('kappa:                   ', kappa,file=write_metadata)
        print('tau:                     ', tau,file=write_metadata)
        print('number_of_particles:     ', Np**2,file=write_metadata)
        print('frac of active cells:    ', frac,file=write_metadata)
        print('gEE:                     ', gEE,file=write_metadata)
        print('gEM:                     ', gEM,file=write_metadata)
        print('gMM:                     ', gMM,file=write_metadata)
        
        print('\nCell types',file=write_metadata)
        print('Medium:                  type = 0',file=write_metadata)
        print('Epithelial:              type = 1',file=write_metadata)
        print('Mesenchymal:             type = 2',file=write_metadata)
        
        temp = self.get_xml_element('temp')
        JMedE = self.get_xml_element('JMedE')
        JMedM = self.get_xml_element('JMedM')
        JEE = self.get_xml_element('JEE')
        JEM = self.get_xml_element('JEM')
        JMM = self.get_xml_element('JMM')
        neighorder= self.get_xml_element('neighorder')
        max_mcs = self.get_xml_element('numbersteps')
        
        
        print('\nXML parameters',file=write_metadata)
        print('box-size x:              ', self.dim.x,file=write_metadata)
        print('box-size y:              ', self.dim.y,file=write_metadata)
        print('box-size z:              ', self.dim.z,file=write_metadata)
        print('J_MedE:                  ',float(JMedE.cdata),file=write_metadata) 
        print('J_MedM:                  ', float(JMedM.cdata),file=write_metadata) 
        print('J_EE:                    ', float(JEE.cdata),file=write_metadata) 
        print('J_EM:                    ', float(JEM.cdata),file=write_metadata) 
        print('J_MM:                    ', float(JMM.cdata),file=write_metadata) 
        print('neighbor order:          ', float(neighorder.cdata),file=write_metadata) 
        print('total mcs:               ', float(max_mcs.cdata),file=write_metadata) 

        print('temperature:             ', float(temp.cdata),'\n',file=write_metadata) 

                
        print('time stamp:              {}h{}m{}s'.format(localtime()[3],localtime()[4],localtime()[5]),file=write_metadata)        
        print('day stamp:               {}d{}m{}y'.format(localtime()[2],localtime()[1],localtime()[0]),file=write_metadata)   

        

    def step(self,mcs):
        if mcs == 2500:
            for cell in self.cell_list:
                cell.dict['angle_pol_now'] = np.pi*random.uniform(-1,1)
                
                cell.dict['angle_pol_update'] = cell.dict['angle_pol_now']
           
               

        if mcs > 2500:
            for cell in self.cell_list:
                Pix = np.cos(cell.dict['angle_pol_now'])
                Piy = np.sin(cell.dict['angle_pol_now'])
                
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if neighbor.type == self.MESENCHYMAL and cell.type == self.MESENCHYMAL:
                            Pix += gMM*np.cos(neighbor.dict['angle_pol_now'])
                            Piy += gMM*np.sin(neighbor.dict['angle_pol_now'])
                
                theta = np.arctan2(Piy,Pix) + 1 * (np.sqrt(2./(tau))*np.random.normal(loc = 0.0, scale=1.0))
                cell.dict['angle_pol_update'] = theta
        
            for cell in self.cell_list:
                if cell.type == self.MESENCHYMAL:
                    cell.lambdaVecX = -kappa*np.cos(cell.dict['angle_pol_update'])
                    cell.lambdaVecY = -kappa*np.sin(cell.dict['angle_pol_update'])
                    cell.lambdaVecZ = 0
                    cell.dict['angle_pol_now'] = cell.dict['angle_pol_update']




        if mcs>=2500 and mcs%200 == 0:     

            ## creating neighborlists
            cell_and_its_neighbors = []
            mesenchymal_cell_id = []
            for cell in self.cell_list_by_type(self.MESENCHYMAL):
                mesenchymal_cell_id.append(cell.id)
            
                neighbor_list = []
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor and neighbor.type == self.MESENCHYMAL:
                        neighbor_list.append(neighbor.id)
                        
                        
                cell_and_its_neighbors.append([cell.id,neighbor_list])
                
                
            cluster_check_list = []
            cluster_id = 0
            cluster_size = 0
            
            cluster_distribution = []
            cluster_id_per_cell = []
            cell_id_already_checked = []
            
            while len(mesenchymal_cell_id)>0:
                
                ### check if there is any element in the cluster_check list. If not, start new cluster
                if len(cluster_check_list) == 0:
                    
                    cluster_id +=1
                    cellid = mesenchymal_cell_id[0] ## choose a new cell to check the cluster size of the cluster it is part of
                    cluster_size = 1
                    
                    ### check neighbors and fill cluster_check_list if applicable
                    for i in range(len(cell_and_its_neighbors)):
                        if cell_and_its_neighbors[i][0] == cellid:
                            cluster_check_list += cell_and_its_neighbors[i][1] ## add neighbors, i.e. cells in the cluster, to the cluster_list that needs to be checked
                           
                    ## keep track of cells that have been considered already
                    mesenchymal_cell_id.remove(cellid)
                    cell_id_already_checked.append(cellid)       
                    
                    if len(cluster_check_list) == 0:
                        ## if a cluster has no more new neighbors, all cells in the cluster have been included and the cluster is closed
                        cluster_distribution.append([cluster_id, cluster_size])                
                    
                    cluster_id_per_cell.append([cellid,cluster_id])
                    
                else: 
                    cellid = cluster_check_list[0] ## choose a cell that is part of an existing cluster
                    cluster_size += 1
                    
                    ### check neighbors and fill cluster_check_list if applicable
                    for i in range(len(cell_and_its_neighbors)):
                        if cell_and_its_neighbors[i][0] == cellid:
                            cluster_check_list += cell_and_its_neighbors[i][1]

                    ## keep track of cells that have been considered already
                    mesenchymal_cell_id.remove(cellid)
                    cell_id_already_checked.append(cellid)                

                    cluster_check_list = list(set(cluster_check_list))
                    for checked_id in cell_id_already_checked:
                        ## remove all cells that have been checked 
                        if checked_id in cluster_check_list:
                            cluster_check_list.remove(checked_id)
                    
                    if len(cluster_check_list) == 0:
                        ## if a cluster has no more new neighbors, all cells in the cluster have been included and the cluster is closed
                        cluster_distribution.append([cluster_id, cluster_size])  
                    cluster_id_per_cell.append([cellid,cluster_id])
            
            for i in range(len(cluster_distribution)):
                print( '{} {} {} '.format(mcs,cluster_distribution[i][0],cluster_distribution[i][1]),file=write_clusters)
              
            for cell in self.cell_list:
                if cell.type == self.EPITHELIAL:
                    cluster_id_cell = 0
                if cell.type == self.MESENCHYMAL:
                    for i in range(len(cluster_id_per_cell)):
                        if cluster_id_per_cell[i][0] == cell.id:
                                cluster_id_cell = cluster_id_per_cell[i][1]
                print( '{} {} {} {} {} {} {} {}'.format(mcs, cell.id, cell.type, cell.xCOM, cell.yCOM, cell.zCOM, cell.dict['angle_pol_now'], cluster_id_cell),file=write_cells)    

            
    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        write_clusters.close()
        write_cells.close()

    def on_stop(self):
        # this gets called each time user stops simulation
        return