#!/bin/bash

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$scriptdir/../.." && pwd)"
outdir="${project_root}/data/SSC_virtualmoorings"
logfile="${outdir}/download_done.log"

mkdir -p "$outdir"
touch "$logfile"

location_name="$1"
X="$2"
Y="$3"
echo "location_name = [$location_name]"
echo "X = [$X]"
echo "Y = [$Y]"

vars=(dissolved_oxygen total_alkalinity dissolved_inorganic_carbon temperature salinity)

years=({2015..2026})
months=({1..12}) 

download_file () {
    var=$1
    year=$2
    month=$3

    mm=$(printf "%02d" "$month")

    outfile="${outdir}/${location_name}_${var}_${year}_${mm}.nc"
    key="${location_name},${var},${year}-${mm}"

    if grep -q "^${key}$" "$logfile"; then
        echo "SKIP (logged): $key"
        return 0
    fi

    if [ -s "$outfile" ]; then
        echo "SKIP (exists): $outfile"
        echo "$key" >> "$logfile"
        return 0
    fi

    last_day=$(date -j -v+1m -v-1d \
        -f "%Y-%m-%d" \
        "${year}-${mm}-01" \
        +"%d")

    case "$var" in 
	dissolved_oxygen|total_alkalinity|dissolved_inorganic_carbon)
    	    url="https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DChemistryFields1hV21-11.nc?${var}%5B(${year}-${mm}-01T00:30:00Z):1:(${year}-${mm}-${last_day}T23:30:00Z)%5D%5B(0.5000003):1:(441.4661)%5D%5B(${Y}):1:(${Y})%5D%5B(${X}):1:(${X})%5D"
	    ;;
	temperature|salinity)
	    url="https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DPhysicsFields1hV21-11.nc?${var}%5B(${year}-${mm}-01T00:30:00Z):1:(${year}-${mm}-${last_day}T23:30:00Z)%5D%5B(0.5000003):1:(441.4661)%5D%5B(${Y}):1:(${Y})%5D%5B(${X}):1:(${X})%5D"
	    ;;
    esac        

    echo $url

    tmpfile="${outfile}.tmp"

    curl --fail --silent --show-error --retry 6 --retry-delay 10 --retry-all-errors --connect-timeout 60 --max-time 3600 -o "$tmpfile" "$url"

    if [ ! -s "$tmpfile" ]; then
        rm -f "$tmpfile"
        return 1
    fi

    mv "$tmpfile" "$outfile"

    echo "$key" >> "$logfile"
    echo "DONE: $outfile"
}

export -f download_file
export outdir
export logfile
export X
export Y
export location_name

parallel -j 1 download_file ::: "${vars[@]}" ::: "${years[@]}" ::: "${months[@]}"
