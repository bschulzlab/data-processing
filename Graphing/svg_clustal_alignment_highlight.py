import svgwrite as sw
import pandas as pd
import re

# Script for writing an SVG image of a clustal alignment
padding_right = 5  # Setting the padding space on the right of the graph
padding_top = 20  # Setting the padding space on the top of the graph
sequence_start = padding_right + 45  # Setting the position where the sequence would actually be drawn from the right of the graph
max_display_len = 65  # Setting the maximum number of aa display on each line
font_size = 4  # Setting the default fontsize
font_family = "Courier"  # Setting the default font
distance_between_aa = 3  # Setting the distance Between the start of each aa character
# start_segment = 21  #Setting the desired amino acid position of the alignment to be drawn
start_segment = 76
skip_end = 363

# Input for color scheme file which define pattern for color annotation
input_color_scheme = r"C:\Users\localadmin\PycharmProjects\svg_alignment\hmwC_repeats_colouringscheme_DE (1).txt"
# Input for the clustal multiple sequence alignment
input_alignment = r"hmw-C_alignment_MUSCLE_final20190807"
# input_legend = r"C:\Users\localadmin\PycharmProjects\svg_alignment\hmw-A_alignment_legend.txt"
output_svg = r"C:\Users\localadmin\PycharmProjects\svg_alignment\hmwC3.svg"


color_dict = {"quality": {}}  # Setting a color dictionary for each amino acid on the alignment that need to be highlight



def get_color(key, color_dict):
    default_color = "white"
    if not color_dict:
        return default_color
    if key in color_dict:
        return color_dict[key]
    else:
        return default_color


def identify_pattern_positions(pattern, source):
    pn = re.compile(pattern.replace(" ", ""))
    for i in pn.finditer(source):
        print(i)

        if len(i.groups()) > 0:
            for m in range(1, len(i.groups()) + 1):
                yield i.start(m), i.end(m)
        else:
            print(i.start())
            yield i.start(), i.end()

def parse_color_position(color_string):
    color_array = color_string.split(";")
    for color in color_array:
        if color:
            yield color.split(",")

with open(
        input_alignment,
        "rt") as align_file:
    sequences = {"quality": ""}
    names = []
    for i in align_file:
        if not i.startswith(" "):
            i = i.strip()
            if not i.startswith("CLUSTAL") and i:
                seq_title = i[0:22].strip()
                if seq_title not in sequences:
                    sequences[seq_title] = ""
                    names.append(seq_title)

                sequences[seq_title] += i[22:].split("\t")[0]
        else:
            sequences["quality"] += i[22:].strip("\n").strip("\r")
    names.append("quality")
    max_seq_len = len(sequences["quality"])

gap_dict = {}

for i in sequences:
    if i != "quality":
        if i not in gap_dict:
            gap_dict[i] = {}
        gap_count = 0
        for p, aa in enumerate(sequences[i]):
            if aa == "-":
                gap_count += 1
            else:
                gap_dict[i][p - gap_count] = gap_count


color_scheme = pd.read_csv(input_color_scheme, sep="\t")
start_dict = {"quality": set()}
end_dict = {"quality": set()}
for i, r in color_scheme.iterrows():
    seq_n = r["Strain"].replace(" ", "_")
    if seq_n not in color_dict:
        color_dict[seq_n] = {}
        start_dict[seq_n] = set()
        end_dict[seq_n] = set()
    for start, end in identify_pattern_positions(r["Sequence"], sequences[seq_n].replace("-", "")):
        start_dict[seq_n].add(start + gap_dict[seq_n][start])
        if ";" in r["Color"]:
            color_array = r["Color"].split(";")
            for color_start, color_end, color_value in parse_color_position(r["Color"]):
                # start_dict[seq_n].add(int(color_start) - 1 + start + gap_dict[seq_n][int(color_start) - 1 + start])
                # end_dict[seq_n].add(int(color_end) + start + gap_dict[seq_n][int(color_end) + start])
                for pos in range(int(color_start) - 1 + start + gap_dict[seq_n][int(color_start) - 1 + start],
                                 int(color_end) + start + gap_dict[seq_n][int(color_end) + start]):

                    if sequences[seq_n][pos] != "-":

                        color_dict[seq_n][pos] = color_value
        else:
            # start_dict[seq_n].add(start + gap_dict[seq_n][start])
            # end_dict[seq_n].add(end + gap_dict[seq_n][end])
            print(start, end)
            for pos in range(start + gap_dict[seq_n][start], end + gap_dict[seq_n][end-1]):

                if sequences[seq_n][pos] != "-":
                    print(pos)
                    color_dict[seq_n][pos] = r["Color"]

# legend_df = pd.read_csv(input_legend, sep="\t")
def draw_aa(group, canvas, aa, pos_x, pos_y, block_size, annotation_color, font_size=4, font_style="normal", font_family="Courier"):
    aa_group = group.add(canvas.g(font_style=font_style, font_family=font_family, font_size=font_size))
    aa_group.add(canvas.rect(insert=(pos_x+0.25, pos_y - 3 - 4.5), size=(block_size[0]-0.75, block_size[1]), stroke=annotation_color, fill=annotation_color))
    aa_group.add(canvas.text(aa, insert=(pos_x, pos_y - 5), fill='black'))


def draw_divider(group, canvas, pos_x, pos_y, size=0.5, color="white"):
    group.add(canvas.line(start=(pos_x - 0.25, pos_y - 8),
                                   end=(pos_x - 0.25, pos_y + 2.5),
                                   stroke_width=size, stroke=color))

def draw_seq_name(group, canvas, name, font_size=4, font_style="normal", font_family="Courier"):
    name_group = group.add(canvas.g(font_size=font_size, font_style=font_style, font_family=font_family))
    name_group.add(canvas.text(name, insert=(padding_right, distance_from_top-5), fill='black'))


print(color_dict)
if __name__ == "__main__":
    canvas = sw.Drawing(output_svg, size=(210*sw.mm, 297*sw.mm))

    seq_name_group = canvas.add(canvas.g())
    distance_from_top = padding_top + 40

    for seg in range(start_segment, max_seq_len, max_display_len):
        top_scale = canvas.add(canvas.g())
        if max_seq_len - max_display_len > seg:
            for i in range(0, max_display_len, 65):
                top_scale.add(
                    canvas.text(str(i + seg + 1 - start_segment),
                                insert=(sequence_start + 1 + i * distance_between_aa, distance_from_top - 15),
                                font_size=font_size,
                                font_family=font_family,
                                fill='black'))
            if i <= max_display_len:
                top_scale.add(
                    canvas.text(str(max_display_len + seg  - start_segment),
                                insert=(sequence_start + 1 + (max_display_len - 1) * distance_between_aa,
                                        distance_from_top - 15), font_size=font_size,
                                font_family=font_family,
                                fill='black'))
        else:
            for i in range(seg, max_seq_len, 65):
                top_scale.add(
                    canvas.text(str(i + 1 - start_segment),
                                insert=(sequence_start + 1 + (i - seg) * distance_between_aa, distance_from_top - 15),
                                font_size=font_size,
                                font_family=font_family,
                                fill='black'))
            if i <= max_seq_len:
                top_scale.add(
                    canvas.text(str(max_seq_len - start_segment),
                                insert=(sequence_start + 1 + (max_seq_len - seg - 1) * distance_between_aa,
                                        distance_from_top - 15),
                                font_size=font_size, font_family=font_family,
                                fill='black'))
        name_group = canvas.add(canvas.g())
        for n in names:
            if n != "quality":
                draw_seq_name(name_group, canvas, n, font_style="normal", font_size=font_size)
            if max_seq_len - max_display_len > seg:
                for i in range(max_display_len):
                    if sequences[n][i]:
                        seq_name_sub_group = seq_name_group.add(canvas.g())
                        draw_aa(seq_name_sub_group,
                                canvas,
                                sequences[n][i+seg],
                                sequence_start + i * distance_between_aa,
                                distance_from_top,
                                (distance_between_aa, distance_between_aa-0.5),
                                get_color(i + seg, color_dict[n]), font_size=font_size)

                        if i + seg in start_dict[n]:
                            draw_divider(seq_name_sub_group, canvas, sequence_start + i * distance_between_aa, distance_from_top -5)

            else:
                for i in range(seg, max_seq_len):
                    if sequences[n][i]:
                        seq_name_sub_group = seq_name_group.add(canvas.g())
                        draw_aa(seq_name_sub_group,
                                canvas,
                                sequences[n][i],
                                sequence_start + (i - seg) * distance_between_aa,
                                distance_from_top,
                                (distance_between_aa, distance_between_aa -0.5),
                                get_color(i, color_dict[n]), font_size=font_size)

                        if i in start_dict[n]:
                            draw_divider(seq_name_sub_group, canvas, sequence_start + (i - seg) * distance_between_aa,
                                         distance_from_top-5)

            distance_from_top += 7
        distance_from_top += 15

    # legend = canvas.add(canvas.g(id="legend", class_="legend"))
    # for i, r in legend_df.iterrows():
    #     legend_color_dict = {}
    #
    #     if ";" in r["Color"]:
    #         for color_start, color_end, color_value in parse_color_position(r["Color"]):
    #             for pos in range(int(color_start)-1, int(color_end)):
    #                 legend_color_dict[pos] = color_value
    #     else:
    #         if not pd.isnull(r["Sequence"]):
    #             legend_color_dict = {k: r["Color"] for k in range(len(r["Sequence"]))}
    #         else:
    #             legend_color_dict = {0: r["Color"]}
    #     if not pd.isnull(r["Sequence"]):
    #         for aa in range(len(r["Sequence"])):
    #             legend.add(canvas.rect(insert=(sequence_start + aa * distance_between_aa - 2, distance_from_top - 8.5),
    #                         size=(10, 10), stroke=legend_color_dict[aa],
    #                         fill=legend_color_dict[aa]))
    #             legend.add(canvas.text(r["Sequence"][aa], insert=(
    #                         sequence_start + aa * distance_between_aa, distance_from_top), font_size=font_size,
    #                                                    font_family=font_family,
    #                                                    fill='black'))
    #     else:
    #         legend.add(canvas.rect(insert=(sequence_start + 0 * distance_between_aa - 2, distance_from_top - 8.5),
    #                                size=(10, 10), stroke=legend_color_dict[0],
    #                                fill=legend_color_dict[0]))
    #         legend.add(canvas.text(r["Description"], insert=(
    #             sequence_start + distance_between_aa, distance_from_top), font_size=font_size,
    #                                font_family=font_family,
    #                                fill='black'))
    #     distance_from_top += 20
    canvas.save()
