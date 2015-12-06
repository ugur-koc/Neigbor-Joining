DATA=85VASTdomains

INPUT="data/${DATA}_distances.txt"
OUTPUT="data/${DATA}_tree.txt"
echo $INPUT
echo $OUTPUT
./main $INPUT $OUTPUT
