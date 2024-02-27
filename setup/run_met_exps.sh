#run burn and unburned cases in serial
./crproj YKD-burned-5yr && ./launch YKD-burned-5yr 
./crproj YKD-burned-no-snow-5yr && ./launch YKD-burned-no-snow-5yr
./crproj YKD-burned-swin20pct-5yr && ./launch YKD-burned-swin20pct-5yr
./crproj YKD-burned-tplus1pct-5yr && ./launch YKD-burned-tplus1pct-5yr
#unburned case
./crproj YKD-unburned-5yr && ./launch YKD-unburned-5yr
./crproj YKD-unburned-no-snow-5yr && ./launch YKD-unburned-no-snow-5yr 
./crproj YKD-unburned-swin20pct-5yr && ./launch YKD-unburned-swin20pct-5yr
./crproj YKD-unburned-tplus1pct-5yr && ./launch YKD-unburned-tplus1pct-5yr


