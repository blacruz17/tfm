# activate conda environment:
conda activate flan 
# format lefse input:
format_input.py lefse.txt lefse.in -c 2 -u 1 -o 1000000
# run lefse with the input file:
# uses 3.5 as LDA cutoff
run_lefse.py lefse.in lefse.res -l 3.5

# plot results:
# barplot:
plot_res.py lefse.res lefse.bar.png --format png --dpi 300 --title '' --max_feature_len 100 --right_space 0.05 --left_space 0.45 --feature_font_size 8 --width 6.25
# cladogram:
plot_cladogram.py lefse.res lefse.cladogram.png --format png --dpi 300 --title '' --label_font_size 5.5 --expand_void_lev 0 --colored_labels 0 --left_space_prop 0.1 --right_space_prop 0.3 --labeled_stop_lev 5 --abrv_start_lev 0 --title_font_size 0 --class_legend_font_size 0