import svgwrite as sw
# Require the python package svgwrite
# Input sequence here for graphic generation
sequence = "TTLTNLTTLESIK"
output_filename = "TTLTNLTTLESIK.svg"

def draw_aa(group, canvas, aa, x_pos, y_pos, font_size=10, font_style="normal", font_family="Courier"):
    # Function for drawing amino acid letter
    aa_group = group.add(canvas.g(fill="black", font_style=font_style, font_family=font_family))
    aa_group.add(canvas.text(aa, insert=(x_pos, y_pos), font_size=font_size))


def draw_divider(group, canvas, points, fill="none", stroke="black"):
    # Function for drawing transition divider.
    group.add(canvas.polyline(points, fill=fill, stroke=stroke))


def draw_transition_label(group, canvas, transition, number, pos_x, pos_y, font_size=6, subscript_font_size=3,
                          font_style="normal", font_family="Helvetica"):
    # Function for drawing transition label
    label_group = group.add(canvas.g(fill="black",
                                     font_style=font_style,
                                     font_size=font_size, font_family=font_family))
    transition = label_group.add(canvas.text(transition, insert=(pos_x, pos_y)))
    transition.add(canvas.tspan(number, baseline_shift="sub", font_size=subscript_font_size))


def draw_peptide_transition(group, canvas, sequence, pos_x, pos_y, distance_between_aa=12):
    # Function for drawing the sequence
    sequence_len = len(sequence)
    for i in range(sequence_len):
        # Iterate through each amino acid and draw out
        draw_aa(group, canvas, sequence[i], pos_x + i * distance_between_aa, pos_y)
        if i != sequence_len - 1:
            draw_divider(group, canvas, [(pos_x + i * distance_between_aa + distance_between_aa + 6, pos_y - 8),
                                         (pos_x + i * distance_between_aa + distance_between_aa - 3, pos_y - 8),
                                         (pos_x + i * distance_between_aa + distance_between_aa - 3, pos_y - 3),
                                         (pos_x + i * distance_between_aa + distance_between_aa - 3, pos_y + 2),
                                         (pos_x + i * distance_between_aa, pos_y + 2)])
            draw_transition_label(group, canvas, "b", i + 1, pos_x + i * distance_between_aa, pos_y + 8)
            draw_transition_label(group, canvas, "y", i + 1, (sequence_len - i) * distance_between_aa, pos_y - 11)


if __name__ == "__main__":
    # Drawing the initial canvas and create the svg file
    canvas = sw.Drawing(output_filename, size=(200, 100))
    # Create the initial graphic object group for holding all the svg elements of the graph.
    peptide_group = canvas.add(canvas.g(id="peptide-transition"))
    # Draw the graph
    draw_peptide_transition(peptide_group, canvas, sequence, 10, 20)

    canvas.save()
