

########################################################
################    Function plot    ###################
########################################################
function_plot () {
fig=results_$1.ps
file=$1_results.txt
# number of burnin models
burnin=$2
# maximum number of iteration for the plot
itmax=$3
# height of the scarp
hs=$4
# max age
agemax=$5

echo ""
echo "---"
echo "Loading file: $file"

echo "Burnin period for histograms: $burnin models"

color=black
if [ -s "$file" ]
then
# Histogram of the number of events
nevmax=15;
gawk -F',' -v burnin=$burnin '{if(NR >= burnin) print $1}' $file | gmt pshistogram  -JX10/10 -R0/$nevmax/0/100 -W1,$color -Ba2f1:'N event':/a20f10:'Frequency (\%)':neSW  -Y48 -Z1 -L -K >  $fig
# Number of events over iterations
gawk -F',' '{print NR,$1}' $file | gmt psxy  -JX10/10 -R0/$itmax/0/$nevmax -W1 -Ba5000f1000:'Accepted models':/a1:'Nb events':neSW  -X13 -W0.1,$color -K -O >>  $fig


# Peri-glacial slip-rate over iterations
SRmax=5;
gawk -F',' '{print NR,$2}' $file | gmt psxy  -JX10/5 -R0/$itmax/0/$SRmax -W1 -Ba5000f1000:'Accepted models':/a1f0.5:'SR (mm/yr)':neSW -Y-8 -W0.1,$color -K -O >>  $fig
# Histogram of peri-glacial slip-rate
HISTOmax=20;
gawk -v burnin=$burnin '{if(NR >= burnin) print $2}' $file | gmt pshistogram  -JX10/5 -R0/$SRmax/0/$HISTOmax -W0.2,$color -Ba1f0.5:'SR (mm/yr)':/a10f5:'Frequency (\%)':neSW -X-13  -Z1 -L -K -O >>  $fig



# RMSw over iteration
RMSmax=100
gmt psxy  -JX10/8 -R0/$itmax/0/$RMSmax -L -Ggray -Ba5000f1000:'Accepted models':/a20f10:'RMSw':neSW -X13 -Y-10 -K -O  << EOF >> $fig
0 0
$burnin 0
$burnin 100
0 100
EOF
gawk -F',' '{print NR,$5}' $file | gmt psxy  -J -R -W1,$color  -W0.1,$color -K -O >>  $fig

# Histogram of the RMSw
RMSmax=50
HISTOmin=0
HISTOmax=100
gawk -v burnin=$burnin -F',' '{if(NR >= burnin) print $5}' $file | gmt pshistogram  -JX10/8 -R0/$RMSmax/$HISTOmin/$HISTOmax -W1,$color -Ba10f2:'RMSw':/a20f10:'Frequency (\%)':neSW -X-13 -A -Z1 -L -K -O >>  $fig



# Event age over iterations
gawk  -F',' '{print NR $1 $3}' $file | gawk '{
for (i = 1; i <= $2; ++i) {
col = i+2;
print $1,$col/1000;
}
}'  | gmt psxy  -JX18/15 -R0/$itmax/0/$agemax -Ba5000f1000:'Accepted models':/a2f1:'Age':neSW  -Sc0.01 -G$color -X26 -Y13 -K -O >>  $fig

# Histogram of the inferred AGE
HISTOmin=0
HISTOmax=10
gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print $1 $3}}' $file | gawk   '{
for (i = 1; i <= $1; ++i) {
col = i+1;
if($col>0) print $col/1000;
}
}' | gmt pshistogram  -JX5/15 -R0/$agemax/$HISTOmin/$HISTOmax -W0.2,$color -Ba2f0.5:'Ages ':/a5f2:'Frequency (\%)':nESw -A -X19 -Z1 -L -K -O >>  $fig



# Slip over iterations
gawk -F',' '{print NR $1 $4}' $file | gawk '{
sum = 0;
for (i = 1; i <$2; ++i) {
col = i+2;
sum += $col;
print $1, sum/100;
}
#print ">";
}'  | gmt psxy  -JX18/15 -R0/$itmax/0/$hs -Ba5000f1000:'Accepted models':/a2f1:'Height on the scarp (m)':neSW  -Sc0.01 -G$color -X-19 -Y-19 -K -O >>  $fig

# Histogram of the inferred SLIP
HISTOmin=0
HISTOmax=10

gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print NR $1 $4}}' $file | gawk '{
sum = 0;
for (i = 1; i <$2; ++i) {
col = i+2;
sum += $col;
print sum/100;
}
#print ">";
}'  | gmt pshistogram  -JX5/15 -R0/$hs/$HISTOmin/$HISTOmax -W0.2,$color -Ba2f1:'Height on the scarp (m) ':/a5f1:'Frequency (\%)':nESw -A -X19 -Z1 -L -K -O >>  $fig

## SLIP OVER THE TIME PLOT ##
# Slip over the time
gawk -v burnin=$burnin -F',' '{if(NR>=burnin)print $1 $3 $4}' $file | gawk '{
slip_sum = 0;
print 0, 0;
for (i = 1; i <=$1; ++i) {

col_age = i+1;
col_slip = i+1+$1;
slip = $col_slip;
age = $col_age;

print age/1000, slip_sum/100;
slip_sum += slip;
print age/1000, slip_sum/100;
}
print ">";
}'  | gmt psxy  -JX18/15 -R0/$agemax/0/$hs -Ba5f1:'Age (kyr)':/a2f1:'Cumulative slip (m)':neSW  -W0.1,$color@95 -X-19 -Y-17 -K -O >>  $fig

# Plot the best model
imin=$(gawk -F',' '
BEGIN { minmis=1e10; }
{ if($5<minmis){minmis=$5;imis=NR} }
END   { print imis }'  $file)

minmis=$(gawk -F',' '
BEGIN { minmis=1e10; }
{ if($5<minmis){minmis=$5;imis=NR} }
END   { print minmis }'  $file)
echo "Best model: model number = $imin,  misfit = $minmis"

gawk -v imin=$imin -F',' '{if(NR==imin)print $1 $3 $4}' $file | gawk '{
slip_sum = 0;
print 0, 0;
for (i = 1; i <=$1; ++i) {

col_age = i+1;
col_slip = i+1+$1;
slip = $col_slip;
age = $col_age;

print age/1000, slip_sum/100;
slip_sum += slip;
print age/1000, slip_sum/100;
}
print ">";
}'  | gmt psxy  -JX18/15 -R0/$agemax/0/$hs -Ba5f1:'Age (kyr)':/a2f1:'Cumulative slip (m)':neSW  -W0.1,green  -K -O >>  $fig
gmt pstext -R -J -N -K -O  << EOF >> $fig
15 2  Best model: it=$imin, misfit=$minmis
EOF

# plot the histogram of cumulative slip
HISTOmax=20

gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print NR $1 $4}}' $file | gawk '{
sum = 0;
for (i = 1; i <$2; ++i) {
col = i+2;
sum += $col;
print sum/100;
}
#print ">";
}'  | gmt pshistogram  -JX5/15 -R0/$hs/0/$HISTOmax -W0.2 -G$color -Ba2f1:'':/a5f1:'Frequency (\%)':nESw -A -X19  -Z1 -L -K -O >>  $fig

# plot the histogram of event ages
HISTOmax=10

gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print $1 $3}}' $file | gawk   '{
for (i = 1; i <= $1; ++i) {
col = i+1;
if($col>0) print $col/1000;
}
}' | gmt pshistogram  -JX18/4 -R0/$agemax/0/$HISTOmax -W0.2 -G$color -Ba5f1:'':/a5f2:'Frequency (\%)':nESw  -X-19 -Y-6 -Z1 -L -O >>  $fig


gmt psconvert $fig -Au -P -Tg -E300
else
echo " file is empty "
fi

echo "Plots are located in $fig"
} # end of function plot




########################################################
################    Merging Function ###################
########################################################
# function to merge the results files
# Copy the content of the given file at the end of the 
#           "merged_results.txt"
# input: number of the inversion, number of models of 
#           the burnin period
function_merge () {
file=$1_results.txt
burnin=$2
echo "-- adding chain $1: -> burnin period: $burnin models"

if [ -s "$file" ]
then
#append the merged results file
gawk -v burnin=$burnin '{if(NR>=burnin){print $0 >> "merged_results.txt" }}' $file
fi
} # end of merging function



########################################################
################    Function grid    ###################
########################################################
# Function to plot the 2D grid density function 
#       (2D frequency histogram)
function_grid () {
fig=grid_$1.ps
file=$1_results.txt

rm frequency.txt -f
rm frequency_5.txt -f
rm frequency_95.txt -f
rm frequency_50.txt -f
rm frequency_ratio.txt -f
rm frequency_sorted.txt -f

burnin=0;
hmax=$2;
Tmax=$3;
echo "Burnin period for histograms: $burnin models"
itmax=50000
color=black
if [ -s "$file" ]
then
# Slip over the time grid
rm -f tmp.txt
gawk -v burnin=$burnin -F',' '{if(NR>=burnin)print $1 $3 $4}' $file > tmp.txt
echo "Gridding"

gawk -v hmax=$hmax -v Tmax=$Tmax '
BEGIN {
Z=0;
dz=10;
T=0;
dt=100;


#INITIZALITION OF THE HISTORGRAM

while(Z<=hmax){
k=k+1;
T = 0;
m=0;
while(T<=Tmax){
m=m+1;
hist[m,k]=0; 
T=T+dt;}
Z=Z+dz;} # EOF while loop
} # EOF BEGIN

{
Z=0;
T = 0;

slip_sum = 0;
age_bef = 0
m=0;
k=0;
for (i = 1; i <=$1; ++i) {
col_age = i+1;
col_slip = i+1+$1;
slip = $col_slip;
age = $col_age;

slip_bef = slip_sum;
slip_sum += slip;
slip_aft = slip_sum;

while(T<=age){
m=m+1;
hist[m,k]=hist[m,k]+1;
#print "-T:",T,"Z",Z,"m:",m,"k:",k,"hist:", hist[m,k];
T=T+dt;
}

while(Z<=slip_aft){
k=k+1;
hist[m,k]=hist[m,k]+1;
#print "T:",T,"-Z",Z,"m:",m,"k:",k,"hist:", hist[m,k];
Z=Z+dz;} # EOF while loop

} # EOF for loop



}

END{
Z=0;
T=0;
k=0;
m=0;
while(Z<=hmax){
k=k+1;
T = 0;
m=0;
Tsum[k]=0;
while(T<=Tmax){
m=m+1;
Tsum[k]=Tsum[k]+hist[m,k];
T=T+dt;}
Z=Z+dz;} # EOF while loop

#writing frequency histogram
Z=0;
T=0;
k=0;
m=0;
hist_cum=0;
while(Z<=hmax){
k=k+1;
T = 0;
m=0;
histrel_max=0;
while(T<=Tmax){
m=m+1;
#            histrel=hist[m,k]/Tsum[k]*100; # relative histogram with frequency in percent
histrel=hist[m,k]; # relative histogram with frequency in percent
hist_cum=hist_cum+hist[m,k];
print T,Z,histrel >> "frequency.txt";
T=T+dt;}
Z=Z+dz;} # EOF while loop
}' tmp.txt

#Sorting histogram by descending frequency
sort -r -nk3 frequency.txt > frequency_sorted.txt

# number of models
n_models=$(gawk '
BEGIN{
nb_model=0;
}
{nb_model=nb_model+$3;
}
END{
print nb_model;
}' frequency_sorted.txt)

# select the n more frequent models in the histogram
ratio=0.95
n_more_freq=$(echo "$n_models*$ratio" | bc)

gawk -v nmax=$n_models '
BEGIN{nb_model=0;}
{nb_model=nb_model+$3;
ratio=(nb_model/nmax)*100
print $1,$2,ratio >> "frequency_ratio.txt";
}' frequency_sorted.txt


gawk '{print $1,$2,$3}' frequency.txt | gmt xyz2grd -R0/$Tmax/0/$hmax -I100/10  -Ghisto.nc
zmax=$(gmt grdinfo -T histo.nc | gawk  -F'/' '{print $2}')
gmt grdmath histo.nc $zmax DIV 100 MUL = histo_rel.nc
gmt grdmath histo.nc 10000 DIV = histo2.nc

gawk '{print $1,$2,$3}' frequency_ratio.txt | gmt xyz2grd -R0/$Tmax/0/$hmax -I100/10  -Gratio.nc

gmt gmtset COLOR_NAN=white
gmt gmtset COLOR_BACKGROUND=white
gmt gmtset COLOR_FOREGROUND=white
gmt grd2cpt histo.nc  -Cgray -M -I -Z > histo.cpt
gmt grd2cpt histo2.nc  -Cgray -M -I -Z > histo2.cpt

gmt grdimage histo.nc -JX20/10 -R0/$Tmax/0/$hmax -Ba5000f1000/a200f20:'Cumulative slip (m)':nesW -X15 -Y15 -K -Chisto.cpt  > $fig
gmt psscale -F+gwhite -D16/3/5/1h -Chisto2.cpt -B10:"Density (x10^4)": -O -K -P >> $fig
gmt grdcontour histo_rel.nc -J -C+5 -S100 -W0.01p,red -A- -O -K  >> $fig


#Export Slip histogram
gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print NR $1 $4}}' $file | gawk '{
sum = 0;
for (i = 1; i <$2; ++i) {
col = i+2;
sum += $col;
print sum;}
}'  | gmt pshistogram  -JX3/20 -R0/$hmax/0/180000 -W10 -G$color  -A -Z0 -IO -L -K -F >  histo_slips.txt
#Export Ages histogram
gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print $1 $3}}' $file | gawk   '{
for (i = 1; i <= $1; ++i) {
col = i+1;
if($col>0) print $col;}
}' | gmt pshistogram  -JX20/3 -R0/$Tmax/0/180000 -W200 -G$color -Z0 -IO  -L -O >  histo_ages.txt

#Sorting histogram by descending frequency
sort -r -nk2 histo_ages.txt > histo_ages_sorted.txt
sort -r -nk2 histo_slips.txt > histo_slips_sorted.txt
# number of models
n_models=$(gawk '
BEGIN{
nb_model=0;
}
{nb_model=nb_model+$2;
}
END{
print nb_model;
}' histo_ages_sorted.txt)

# select the n more frequent models in the Age histogram
rm -rf histo_ages_ratio.txt histo_slips_ratio.txt

ratiomax=0.95 # We keep 95% of the most frequent models
n_more_freq=$(echo "$n_models*$ratio" | bc)

# Selecting the 95% most frequent ages
gawk -v nmax=$n_models -v ratiomax=$ratiomax '
BEGIN{nb_model=0;}
{nb_model=nb_model+$2;
ratio=(nb_model/nmax)
if(ratio<=ratiomax){
print $1,$2 >> "histo_ages_ratio.txt";}
}' histo_ages_sorted.txt

# number of models
n_models=$(gawk '
BEGIN{
nb_model=0;
}
{nb_model=nb_model+$2;
}
END{
print nb_model;
}' histo_slips_sorted.txt)


ratiomax=0.95 # We keep 95% of the most frequent models
n_more_freq=$(echo "$n_models*$ratio" | bc)

# Selecting the 95% most frequent slips
gawk -v nmax=$n_models -v ratiomax=$ratiomax '
BEGIN{nb_model=0;}
{nb_model=nb_model+$2;
ratio=(nb_model/nmax)
if(ratio<=ratiomax){
print $1,$2 >> "histo_slips_ratio.txt";}
}' histo_slips_sorted.txt

#Plot Slip histogram
max_density=600000;
max_density2=$(echo "$max_density/10000" | bc);
gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print NR $1 $4}}' $file | gawk '{
sum = 0;
for (i = 1; i <$2; ++i) {
col = i+2;
sum += $col;
print sum;}
}'  | gmt pshistogram  -JX3/10 -R0/$hmax/0/$max_density -W10 -Ggrey  -A -X21  -Z0 -L -K -F -O >>  $fig
#Plot the 95% slips
#gawk '{print $2/10000, $1+5}' histo_slips_ratio.txt | gmt psxy  -J -R0/$max_density2/0/$hmax -Ba20f10:'Density (x10^4)':/a200f100:'':nESw -SB10u -Gred@50 -K  -O  >>  $fig

#Plot Ages histogram
max_density=400000;
max_density2=$(echo "$max_density/10000" | bc)
gawk -v burnin=$burnin -F',' '{if(NR>=burnin){print $1 $3}}' $file | gawk   '{
for (i = 1; i <= $1; ++i) {
col = i+1;
if($col>0) print $col;}
}' | gmt pshistogram  -JX20/3 -R0/$Tmax/0/$max_density -W200 -Ggrey  -X-21 -Y-4 -Z0 -K -L -O >>  $fig
#Plot the 95% ages
#gawk '{print $1+100, $2/10000}' histo_ages_ratio.txt | gmt psxy  -J -R0/$Tmax/0/$max_density2 -Ba5000f1000:'Age (kyr)':/a20f10:'Density (x10^4)':neSW  -Sb200u -Gred@50 -O >>  $fig

gmt psconvert $fig -Au -P -Tg -E300

rm -f density.nc tmp.txt frequency_ratio.txt frequency_sorted.txt

else
echo " file is empty "
fi


rm -f density.nc tmp.txt frequency_sorted.txt frequency.txt frequency_ratio.txt
rm -f histo.nc
rm -f ratio.nc
rm -f histo_rel.nc
rm -f histo2.nc
rm -f histo.cpt
rm -f histo2.cpt
rm -f histo_slips.txt
rm -f histo_ages.txt
rm -f histo_ages_sorted.txt
rm -f histo_slips_sorted.txt
rm -f histo_ages_ratio.txt
rm -f histo_slips_ratio.txt

} # end of function plot

