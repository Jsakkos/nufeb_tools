from pathlib import Path
import json
from datafed.CommandLib import API
from nufeb_tools import __version__

# TODO Fix datafed uploading
def upload_datafed():
    """Collect NUFEB simulation data and upload to a DataFed collection
    """
    metadata_file = Path('metadata.json')
    if metadata_file.is_file():
        f = open(metadata_file)
        metadata  = json.load(f)
        f.close()
        n_cyanos = metadata['cyano']['StartingCells']
        n_ecw = metadata['ecw']['StartingCells']
        dims = metadata['Dimensions']
        SucPct = int(metadata['SucRatio']*100)
    else:
        print('No metadata file')
        metadata=None

    try:
        df_api = API()
        df_api.setContext('p/eng107')
        collectionName = f'NUFEB_{n_cyanos}_{n_ecw}_{SucPct}_{dims[0]}_{dims[1]}_{dims[2]}'
        parent_collection = df_api.getAuthUser().split('/')[1]
        coll_msg = df_api.collectionCreate(collectionName,  parent_id=parent_collection)
        global_coll_id = coll_msg[0].coll[0].id
    #_logger.info(global_coll_id)
    except:
        global_coll_id = None
    file_types = ['VTK.tar.gz','trajectory.h5','nufeb.log',f'Results/avg_concentration.csv',f'Results/biomass.csv',f'Results/ntypes.csv']
    titles = ['VTK','Trajectory','Log','Avg_con','Biomass','Ntypes']
    for file,title in zip(file_types,titles):
        fp = Path(file)
        if fp.is_file():
            filename = file
            file_title= title

            rec_msg = df_api.dataCreate(title = file_title,
                                        alias = '',
                                        metadata=json.dumps(metadata),
                                        parent_id=global_coll_id,
                                            )
            rec_id = rec_msg[0].data[0].id
            #Use as pathname the path and name of the file you wish to move from CADES to DataFed
            pput_msg = df_api.dataPut(rec_id, filename, wait=False)
            #_logger.info(pput_msg)
        #else:
            #_logger.debug('No metadata file found')
        #  sys.exit(1)
        else:
            print('Can not find ' + file)
    print('Created collection ' + global_coll_id)
 
    
def run():
    """Calls :func:`upload_datafed` 

    """
    upload_datafed()


if __name__ == "__main__":

    run()