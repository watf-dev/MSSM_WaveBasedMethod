#!/bin/sh
# Created: Dec, 18, 2025 09:53:16 by Wataru Fukuda
set -eu

mkdir -p figs

for i in {1..9}; do
  data=EVAL_NODE_nGP500/eval_node${i}.txt

  read ymin_ ymax_ <<< $(awk 'NR==1 {min=$2; max=$2} { if ($2 < min) min=$2; if ($2 > max) max=$2 } END { print min, max }' $data)

  range=$(awk -v a=$ymax_ -v b=$ymin_ 'BEGIN{print a-b}')
  ymin=$(awk -v y=$ymin_ -v r=$range 'BEGIN{print y-0.1*r}')
  ymax=$(awk -v y=$ymax_ -v r=$range 'BEGIN{print y+0.1*r}')

  echo node $i: ymin=$ymin ymax=$ymax

  lualatex -interaction=batchmode -shell-escape "\def\datanum{$i} \def\ymin{$ymin} \def\ymax{$ymax} \input{node.tex}"
  pdfcrop node.pdf figs/node${i}.pdf --margin 5
done

rm *.aux *.log *.pdf
