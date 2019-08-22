#/bin/ksh
#extract metadata from the xml file
ver=3.5.5
granule=$2

dir=$1/$granule
file=`ls $dir/*xml`

#GEOMETRICAL CONDITIONS
sza=`xml_grep --cond Mean_Sun_Angle/ZENITH_ANGLE $file --text_only `
saa=`xml_grep --cond Mean_Sun_Angle/AZIMUTH_ANGLE $file --text_only `
echo $sza $saa
vza=`xml_grep --cond Mean_Viewing_Incidence_Angle/ZENITH_ANGLE $file --text_only`  
vaa=`xml_grep --cond Mean_Viewing_Incidence_Angle/AZIMUTH_ANGLE $file --text_only`
vza=`echo $vza | awk '{print $1}'`
vaa=`echo $vaa | awk '{print $1}'`
echo $vza $vaa

#UTM ZONE
utmzone=`xml_grep --cond HORIZONTAL_CS_NAME $file --text_only`
utmzone=`echo $utmzone | awk '{print $5}'`
hemi=`echo -n $utmzone | tail -c -1`
if [[ $hemi == "N" ]] 
then
utmzone=`echo $utmzone | tr -d N`
fi
if [[ $hemi == "S" ]] 
then
utmzone=`echo $utmzone | tr -d S`
utmzone=-$utmzone
fi
echo $utmzone
x0=`xml_grep --cond ULX $file --text_only`  
x0=`echo $x0 | awk '{print $1}'`
y0=`xml_grep --cond ULY $file --text_only`  
y0=`echo $y0 | awk '{print $1}'`
echo $x0 $y0

#EXTRACTION DAY AND YEAR
date=`xml_grep --cond SENSING_TIME $file --text_only`
echo $date
year4d=`echo $date | cut -c 1-4`
year2d=`echo $date | cut -c 3-4`
month=`echo $date | cut -c 6-7`
day=`echo $date | cut -c 9-10`
echo $year4d $year2d $month $day
doy=`jdoy/jdoy $month/$day/$year2d`
doy=`echo $doy | cut -c 1-3`

#READING ANCILLARY DATA (WV and O3 from MODIS)
fileanc=`ls /usr/local/ledaps/L8ANC/LADS/$year4d/L8ANC$year4d$doy*.hdf_fused`
echo $fileanc

#CREATING HDF FILE
#for jpfile in `ls $dir/IMG_DATA/*T32TLP*_B08.jp2`
#do
#echo $jpfile
#nfilehdf=`echo $file | awk -F / '{print $NF}' |sed -e "s/.jp2/.hdf/g" | sed -e "s/_B08//g"`
#echo $nfilehdf
#rm -f $nfilehdf
#dir=`echo $jpfile | awk -F / '{printf "%s/%s/%s/%s/\n",$1,$2,$3,$4}'`
#echo $dir
nfilehdf=$dir/$granule.hdf

#CONVERTING jp2 into HDF
if [ -e "$nfilehdf" ] 
then
echo "the file " $nfilehdf " already exists skipping conversion "
else
echo "we must convert to hdf ... "
for band in 01 02 03 04 05 06 07 08 8A 09 10 11 12
do
#gdal_translate -of GTiff  $dir/IMG_DATA/*_B$band.jp2 temp.tif
gdalwarp -overwrite -of GTiff  $dir/IMG_DATA/*_B$band.jp2 $dir/temp.tif
echo  $dir/IMG_DATA/*_B$band.jp2 
tiff2hdf/tiff2hdf $dir/temp.tif $dir/temp.hdf
SDS_transfer_rename/SDS_transfer_rename $dir/temp.hdf $nfilehdf "BAND 0" "toaband $band"
rm $dir/temp.tif
done
fi
#exit

#done

#CREATING INPUT FILE for MSIprocV2.1
fileinp=$dir/${granule}.inp
echo 0 >$fileinp
echo $nfilehdf >>$fileinp
echo $fileanc >>$fileinp
#xml_grep --cond NROWS $file --text_only >${fileinp}.tmp.nwinp1
#xml_grep --cond NCOLS $file --text_only >${fileinp}.tmp.nwinp2
#paste ${fileinp}.tmp.nwinp1 ${fileinp}.tmp.nwinp2 >>$fileinp
nrows=`xml_grep --cond NROWS $file --text_only`
ncols=`xml_grep --cond NCOLS $file --text_only`
for xx in $nrows
do
  echo ${xx/.*} >> ${fileinp}.tmp.nwinp1
done
for xx in $ncols
do
  echo ${xx/.*} >> ${fileinp}.tmp.nwinp2
done
paste ${fileinp}.tmp.nwinp1 ${fileinp}.tmp.nwinp2 >>$fileinp
echo $sza $saa $vza $vaa >>$fileinp
echo $utmzone 1 1 $y0 $x0 10. >>$fileinp
echo $doy >>$fileinp
echo ${dir}/srLaSRCS2AV${ver}-${granule}.hdf >>$fileinp
#echo  1 1000 1 3000 >>$fileinp

#LAUNCH AC CODE
time ./LaSRCS2AV${ver}  <${fileinp} > $dir/srLaSRCS2AV${ver}-$granule.log

rm ${fileinp}.tmp.nwinp1 ${fileinp}.tmp.nwinp2 $dir/temp.hdf
