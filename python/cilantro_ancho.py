#%%
def logger_parse_RLH(input_log):
    """
    This Python script was derived directly from cilantro_ancho_parse.m from Flat Earth, Inc 
        by Von P. Walden, Washington State University

    % Copyright: 2019 Flat Earth Inc
    % Written by: Justin Hadella, Tyler Davis, Trevor Vannoy, Will Tidd

    Here is the documentation from that script:

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FlatEarth - salsa logger - Data Parse Utility
    %
    % This is a utility which parses binary data log files produced by the 
    % salsa logger app for the X1 radar platform.
    %
    % The data log is in the format:
    % [config]
    % [frame_0]
    % [frame_1]
    % [frame_2]
    % ...
    % [frame_N-1]
    % 
    % where...
    % [config]  - 256B of radar configuration data (defined below)
    % [frame_n] - Array of frame data
    %
    % The [config] is a 256-byte section organized as:
    % [MAGIC#]                      - uint32                      
    % [VERSION]                     - uint8
    % [CAPE]                        - uint8
    % [DACMin]                      - int32
    % [DACMax]                      - int32
    % [DACStep]                     - int32
    % [Iterations]                  - int32
    % [PulsesPerStep]               - int32
    % [FrameStitch]                 - int32
    % [ClkDivider]                  - int32
    % [SamplingRate]                - int32
    % [SamplersPerFrame]            - int32
    % [OffsetDistanceFromReference] - float
    % [SampleDelayToReference]      - float
    % [SamplesPerSecond]            - float
    % [SampleDelay]                 - float
    % [...Padding...]               - 194 x uint8
    %
    % TODO -- GPS changes [frame] definitions!
    %
    % The [frame] is organized as:
    % [frame_bytes] - uint32
    % [radar_bytes] - uint32
    % [num_samples] - uint32
    % [Temperature] - float
    % [Radar Data]  - Array of radar counter values (of type uint32_t)
    %
    % Copyright: 2019 Flat Earth Inc
    % Written by: Justin Hadella, Tyler Davis, Trevor Vannoy, Will Tidd
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Args:
        filename (str): filename of the radar log file
    
    Returns:
        gpr (xarray):   xarray object containing the GPR data
    """
    import numpy as np
    from datetime import datetime
    import xarray as xr
    
    def f_DAC2Volt(DAC, config):
        """
        This Python script was derived directly from e_snowdry.m by H.P. Marshall of Boise State University 
        by Von P. Walden, Washington State University

        Here is the documentation from that script:
    
            %This function normalizes the DACcount to a voltage according to the
            %documentation ('Salsa Radar Primer.pdf')
            %DAC = raw data of CSNOW radar
            %config = config file as read in from radar
        
        Args:
            DAC (float): Array of raw data (datar) from the GPR
            config (dict): Python dictionary containing the configuration parameters of the radar

        Returns:
            Volt (float): Array of GPR data as voltages
        """
        normDAC = DAC * config['DACStep']/(config['PulsesPerStep']*config['Iterations'])+config['DACMin']
        Volt = normDAC * 1.04 / 8191
        
        return Volt

    TIMESTAMP_BYTES = 8;
    SWICH_CONFIG_BYTES = 4;

    # ....Opens the DMV file.
    fid = open(input_log, 'rb')

    # ....Determine the file size by searching for the end-of-file; eof.
    eof = fid.seek(-1, 2)  # go to the file end and record byte value
    fileLength = fid.tell()

    # ....Return to beginning of the file/header.
    fid.seek(0)

    # ....Read the file header for configuration parameters
    config = {'magic':                       np.fromfile(fid, np.uint32, 1)[0], 
              'cape':                        np.fromfile(fid, np.uint8,  1)[0],
              'ver':                         np.fromfile(fid, np.uint8,  1)[0],
              'DACMin':                      np.fromfile(fid, np.uint32, 1)[0],
              'DACMax':                      np.fromfile(fid, np.uint32, 1)[0],
              'DACStep':                     np.fromfile(fid, np.uint32, 1)[0],
              'Iterations':                  np.fromfile(fid, np.uint32, 1)[0],
              'PulsesPerStep':               np.fromfile(fid, np.uint32, 1)[0],
              'FrameStitch':                 np.fromfile(fid, np.uint32, 1)[0],
              'PGSelect':                    np.fromfile(fid, np.uint32, 1)[0],
              'SamplingRate':                np.fromfile(fid, np.uint32, 1)[0],
              'SamplersPerFrame':            np.fromfile(fid, np.uint32, 1)[0],
              'OffsetDistanceFromReference': np.fromfile(fid, np.float32,1)[0],
              'SampleDelayToReference':      np.fromfile(fid, np.float32,1)[0],
              'SamplesPerSecond':            np.fromfile(fid, np.float32,1)[0],
              'SampleDelay':                 np.fromfile(fid, np.float32,1)[0],
              'TXVoltage':                   np.fromfile(fid, np.float32,1)[0],
              'NumTrials':                   np.fromfile(fid, np.uint32, 1)[0],
              'FrameSize':                   np.fromfile(fid, np.uint32, 1)[0],
              'NumFrames':                   np.fromfile(fid, np.uint32, 1)[0],
    }
    if config['cape']:
        config['cape'] = 'Ancho'
    else:
        config['cape'] = 'Unknown'
    
    # ....Read (and discard) the rest of the header
    headerJunk = np.fromfile(fid, np.uint8, 182)

    # ....Compute the size of the [frame] MODIFIED RLH TO 1 CH
    config['FrameSize'] = config['SamplersPerFrame']  + TIMESTAMP_BYTES + SWICH_CONFIG_BYTES

    # ....Compute length of radar signal
    signalLength = config['SamplersPerFrame']

    # ....Compute the number of frames in the log (/4 to account for uint32)
    config['NumFrames'] = round((fileLength - 256) / config['FrameSize'] / 4);

    # ....Create default arrays containing NaNs
    tv            = np.nan * np.ones(config['NumFrames'])
    switch_config = np.nan * np.ones(config['NumFrames'])
    datar         = np.nan * np.ones([config['NumFrames'], config['SamplersPerFrame']])

    # ....Read the frames of data
    for scan in range(config['NumFrames']):
        # ....Decode time
        tv_sec  = np.fromfile(fid, np.uint32, 1)[0]
        tv_nsec = np.fromfile(fid, np.uint32, 1)[0]
        tv[scan] = tv_sec + tv_nsec/1e9

        # ....Keep track of switch_config
        switch_config[scan] = np.fromfile(fid, np.uint32, 1)[0]

        # ....Read data
        datar[scan] = np.fromfile(fid, np.int32, config['SamplersPerFrame'])

    # ....Close input_log file
    fid.close()
    
    # ....Create time and step arrays
    time = [datetime.fromtimestamp(t) for t in tv]
    step = np.arange(config['SamplersPerFrame'])

    # ....Convert data from DAC to voltage
    DM = f_DAC2Volt(datar.mean(axis=0), config)
    time_DM = datetime.fromtimestamp(tv.mean())

    # ....Create a netCDF object of GPR data
    gpr = xr.Dataset(
        data_vars=dict(
            datar=(['time', 'step'], datar),
            DM=(['step'], DM),
            time_DM=time_DM
        ),
        coords=dict(
            time=(['time'], time),
            step=(['step'], step)
        ),
        attrs=config,
        )

    return gpr

# %%
def e_snowdry(rho, f=10e9, T=-10):
    """
    This Python script was derived directly from e_snowdry.m by H.P. Marshall of Boise State University 
    by Von P. Walden, Washington State University
    
    Here is the documentation from that script:

        % HPM 09/19/03
        % this function calculates the real and complex part of the dielectric
        %   constant of dry snow, given the density, frequency and temperature
        %   After [Tiuri et al, 1984]
        % INPUT: rho = dry density of snow [kg/m^3]
        %          f = frequency [Hz] (default 10e9)
        %          T = temperature [deg C] (default -10)
        % OUTPUT: e_s = complex diel constant of the snow
    
    Args:
        rho (float):   dry density of snow [kg/m^3]
        f (float):     frequency [Hz] (default 10e9)
        T (float):     temperature [deg C] (default -10)
        
    Returns:
        e_s (complex): complex dielectric constant of dry snow
    """
    import numpy as np
    
    rho=rho/1000                                                # change units to [g/cc]
    e_r=1+1.7*rho+0.7*rho**2                                    # real part of diel const of dry snow
    e_ice=1.59e6*(1./f+1.23e-14*np.sqrt(f))*np.exp(0.036*T)     # imaginary part of diel const of ice, from Tiuri, 1984 eq 6
    e_i=(0.52*rho+0.62*rho**2)*e_ice                            # imaginary part of diel const of dry snow
    e_s=complex(e_r, e_i)                                       # complex diel const of dry snow
    
    return e_s

# %%