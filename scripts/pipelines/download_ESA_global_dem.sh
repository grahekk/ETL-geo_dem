#!/bin/bash
URL="https://prism-dem-open.copernicus.eu/pd-desk-open-access/publicDemURLs/COP-DEM_GLO-90-DGED__2021_1"

# Use curl to fetch the list of download links and then extract them
LINKS=$(curl -s -k -H "accept: xml" "$URL" | grep -oP '(?<=<nativeDemUrl>).*(?=</nativeDemUrl>)')

# Loop through the extracted links and use wget to download each file
for link in $LINKS; do
    echo "link downloaded" $link
    wget -c "$link"
done
