#%% IMPORTING LIBRARIES
import logging
import os
from datetime import datetime

#%% LOGGINF CONFIG
def setup_logging(base_dir='logs', log_dir='optimization_logs'):
    # Create the base directory if it doesn't exist
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    
    # Create the log directory within the base directory
    log_path = os.path.join(base_dir, log_dir)
    os.makedirs(log_path, exist_ok=True)
    
    # Define the log file name with the current timestamp
    log_file = os.path.join(log_path, f'fin_optimization_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logging.info("Logging setup completed.")
    return log_file
# %%
