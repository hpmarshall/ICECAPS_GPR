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
function [data, config, switch_config] = cilantro_ancho_logger_parse_RLH(input_log)
% Outputs: 
% data - matrix containing N radar signals (interleaved according toswitch_config)
% config - configuration structure (to put data into context)
% switch_config - array of the switch settings (4 different settings)
%
% Inputs: input_log - binary data log
%

TIMESTAMP_BYTES = 8;
SWICH_CONFIG_BYTES = 4;

% Open the log file
fid = fopen(input_log);

% Determine total file length
fseek(fid, 0, 'eof');
fileLength = ftell(fid);

% Go back to the beginning
fseek(fid, 0, 'bof');

% The first 256 bytes contain the header which is configuration info
config.magic = fread(fid, 1, 'uint32');
cape = fread(fid, 1, 'uint8');
config.ver = fread(fid, 1, 'uint8');
if cape
  config.cape = 'Ancho';
else
  config.cape = 'Unknown';
end

config.DACMin = fread(fid, 1, 'int32');
config.DACMax = fread(fid, 1, 'int32');
config.DACStep = fread(fid, 1, 'int32');
config.Iterations = fread(fid, 1, 'int32');
config.PulsesPerStep = fread(fid, 1, 'int32');
config.FrameStitch = fread(fid, 1, 'int32');
config.PGSelect = fread(fid, 1, 'int32');
config.SamplingRate = fread(fid, 1, 'int32'); % added this!!
config.SamplersPerFrame = fread(fid, 1, 'int32');
config.OffsetDistanceFromReference = fread(fid, 1, 'float');
config.SampleDelayToReference = fread(fid, 1, 'float');
config.SamplesPerSecond = fread(fid, 1, 'float');
config.SampleDelay = fread(fid, 1, 'float');
config.TXVoltage = fread(fid,1,'float');
config.NumTrials = fread(fid,1,'int32');
fread(fid, 190, 'uint8'); % CHANGED FROM 194 RLH

% Compute the size of the [frame] MODIFIED RLH TO 1 CH
% config.FrameSize = config.SamplersPerFrame * 4 + TIMESTAMP_BYTES + SWICH_CONFIG_BYTES;
config.FrameSize = config.SamplersPerFrame  + TIMESTAMP_BYTES + SWICH_CONFIG_BYTES; 

% Compute length of radar signal
signalLength = config.SamplersPerFrame;

% Compute the number of frames in the log (/4 to account for uint32)
% config.NumFrames = floor((fileLength - 256) / config.FrameSize / 4)
config.NumFrames = round((fileLength - 256) / config.FrameSize / 4);

% Preallocate memory for N radar frames
ts = zeros(1, config.NumFrames);
switch_config = zeros(1, config.NumFrames);
data = zeros(config.NumFrames, signalLength);

for i = 1:config.NumTrials
       % disp(sprintf('%d of %d trials\n',i,config.NumTrials))
       % Read out the timestamp
       tv_sec = fread(fid, 1, 'uint32');
       tv_nsec = fread(fid, 1, 'uint32');
       ts(i) = tv_sec + tv_nsec / 1e9;
       
       % Read out the switch configuration (0,1,2,3,0,1,2,3...)
       switch_config(i) = fread(fid, 1, 'uint32');
       
       % Interleaved data (switch_config 0,1,2,3,0,1,2,3...most likely)
       data(i,:) = fread(fid, [1 signalLength], 'int32');
end

% Close the log file
fclose(fid);






































