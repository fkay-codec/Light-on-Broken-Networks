#!/usr/bin/env bash
# 17 ROI MNI coordinates (x y z), flat array: x1 y1 z1 x2 y2 z2 ...
## we are using this template: MNI152_T1_1mm
channel_id=(
  S1-D1
  S1-D2
  S1-D28
  S2-D1
  S2-D2
  S2-D3
  S2-D4
  S3-D4
  S4-D2
  S4-D3
  S4-D5
  S4-D26
  S5-D3
  S5-D4
  S6-D4
  S7-D3
  S7-D5
  S7-D7
  S8-D4
  S8-D6
  S8-D8
  S9-D5
  S9-D7
  S9-D9
  S9-D23
  S10-D6
  S10-D10
  S11-D8
  S11-D11
  S12-D9
  S12-D10
  S12-D12
  S13-D10
  S14-D9
  S14-D12
  S14-D19
  S15-D10
  S15-D12
  S15-D13
  S15-D14
  S16-D11
  S16-D13
  S17-D12
  S17-D14
  S17-D15
  S17-D18
  S17-D19
  S18-D15
  S18-D17
  S18-D18
  S19-D14
  S19-D15
  S19-D16
  S20-D15
  S21-D18
  S21-D19
  S21-D20
  S21-D21
  S22-D20
  S22-D22
  S23-D9
  S23-D19
  S23-D21
  S24-D21
  S25-D21
  S25-D25
  S26-D22
  S26-D24
  S27-D5
  S27-D23
  S27-D26
  S28-D24
  S28-D25
  S28-D27
  S29-D26
  S29-D27
  S30-D27
  S31-D2
  S31-D26
  S31-D27
  S31-D28
  S32-D27
)

Coords=(
13	67	0
1	64	14
-12	67	0
25	63	9
13	61	24
22	52	33
40	50	16
48	46	5
2	50	39
10	41	50
1	27	58
-9	41	50
30	40	41
46	38	24
55	36	5
24	26	55
14	13	66
27	-4	68
58	24	18
56	12	33
64	-5	22
-1	-4	72
17	-21	75
1	-35	75
-17	-20	74
52	-4	48
53	-35	52
67	-19	4
68	-31	-12
17	-50	73
39	-49	60
25	-62	63
58	-48	38
2	-61	66
15	-73	57
-13	-73	56
46	-62	47
33	-74	48
47	-72	30
34	-83	30
62	-52	-8
54	-67	6
12	-81	43
16	-92	24
-2	-97	12
-16	-93	25
-12	-83	42
-14	-101	-2
-24	-95	-13
-26	-96	8
26	-97	8
15	-99	-1
25	-94	-14
-2	-96	-11
-34	-83	30
-32	-73	47
-46	-72	30
-46	-61	46
-54	-67	6
-64	-52	-9
-16	-50	72
-24	-62	62
-39	-48	60
-57	-48	38
-52	-34	52
-50	-3	50
-68	-32	-12
-65	-18	4
-13	12	67
-26	-5	68
-23	26	56
-62	-3	23
-55	12	34
-56	24	20
-31	39	41
-46	39	26
-53	37	6
-12	62	23
-23	52	32
-39	50	17
-24	63	9
-47	46	6
)

directory="/home/mybox/Desktop/Analysis/OutPuts" # output save
template="/home/mybox/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz"   # adjust path if needed
radius=5	              # sphere radius for -kernel sphere

# Add this before the fslmaths commands:
if [ ! -f "$template" ]; then
    echo "ERROR: Template file not found: $template"
    exit 1
fi

for ((i=0; i<${#Coords[@]}; i+=3)); do
  roi=$(( i/3 + 1 ))
  x_mni=${Coords[i]}
  y_mni=${Coords[i+1]}
  z_mni=${Coords[i+2]}

  channel_index=$((i/3))
  channel_name=${channel_id[channel_index]}
  echo "Channel: $channel_name (index $channel_index)"

  echo "Processing ROI${roi} at MNI coordinates: [$x_mni, $y_mni, $z_mni]"
  # exit

  # MNI -> voxel indices: simple addition
  xi=$((90 - x_mni))  # the correct conversion is VOXEL = 90 - MNI
  yi=$((y_mni + 126)) # VOXEL = MNI + 126
  zi=$((z_mni + 72))  # VOXEL = MNI + 72
  echo "converted to voxel indices: [$xi, $yi, $zi]"
  # exit

  echo "Creating sphere ROI${roi} at voxel [$xi,$yi,$zi] from MNI [$x_mni,$y_mni,$z_mni]"
  # example command: fslmaths avg152T1.nii.gz -mul 0 -add 1 -roi 45 1 74 1 51 1 0 1 ACCpoint -odt float
  echo "$template" -mul 0 -add 1 -roi "$xi" 1 "$yi" 1 "$zi" 1 0 1 "$directory/ROI${roi}_${channel_name}_point" -odt float
  fslmaths "$template" -mul 0 -add 1 -roi "$xi" 1 "$yi" 1 "$zi" 1 0 1 "$directory/ROI${roi}_${channel_name}_point" -odt float # The points are created correctly
  
  # eg. fslmaths ACCpoint -kernel sphere 5 -fmean ACCsphere -odt float
  echo "$directory/ROI${roi}_${channel_name}_point" -kernel sphere "$radius" -fmean "$directory/ROI${roi}_${channel_name}_sphere" -odt float
  fslmaths "$directory/ROI${roi}_${channel_name}_point" -kernel sphere "$radius" -fmean "$directory/ROI${roi}_${channel_name}_${radius}mmSphere" -odt float
  #exit
done 
