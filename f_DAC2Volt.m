function [Volt] = f_DAC2Volt(DAC,config)
%This function normalizes the DACcount to a voltage according to the
%documentation ('Salsa Radar Primer.pdf')
%DAC = raw data of CSNOW radar
%config = config file as read in from radar
normDAC = DAC * config.DACStep/(config.PulsesPerStep*config.Iterations)+config.DACMin;
Volt = normDAC * 1.04 / 8191;
end

