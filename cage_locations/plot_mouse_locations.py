# Plots the locations in protein-carbohydrate intake space of cages in the SSB study. 
# All cages are plotted on the same graph. 
# Each individual's cage is also highlighted in a series of separate graphs. 
#
# Mark N. Read, 2019
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import pandas
from scipy.spatial import ConvexHull
import seaborn as sns


draw_labels = False

# Read SSB input file
meta_fp = '../diet-nutrient-breakdown-SSB.xlsx'
meta_df = pandas.read_excel(io=meta_fp, sheet_name='compound', skiprows=2)

# Subset only the cages of the original SSB study. 
ssb_indices = meta_df['Experiment'] == 'SSB'
meta_ssb_df = meta_df.loc[ssb_indices, :] 
print(meta_ssb_df.shape)  # Sanity check, should be 250. It is. 

# Calculate convex hull of points, to plot. 
# NDArray representation, needed for finding convex hull.
meta_nd = meta_ssb_df.loc[:, ['prot-intake-KJ/d', 'carb-intake-KJ/d']].values
hull = ConvexHull(meta_nd)  # Used below. 

# Colour palette for the diet codes
diet_codes = meta_ssb_df['%P/%C/%F'].unique()
palette = dict(zip(diet_codes, sns.color_palette()))  # Map diet code to palette color. 
for key in palette.keys():
	palette[key] = '#666666'

palette_black = {}  # Using this palette will set all data points to the same colour. 
for key in palette.keys():
    palette_black[key] = 'k'


# All cages together. 
if False:
    ax = sns.scatterplot(x='prot-intake-KJ/d', y='carb-intake-KJ/d', data=meta_ssb_df, 
                        hue='%P/%C/%F', palette=palette, 
                        legend=False, 
                        # linewidth=0,  # Turns off the white marker outline. 
                        s=150, 
                        # style='Energy',  # Gives a different marker shape to each energy level. LOOKS MESSY. 
                        marker='o',  # marker='$\circ$', # Draws a circle, but it's not centred at the coordinate!
                        alpha=0.7)
    # Plot the convex hull of the data.
    for simplex in hull.simplices:
        plt.plot(meta_nd[simplex, 0], meta_nd[simplex, 1], 'k-')

    if draw_labels:
        plt.xlabel('Protein intake (KJ/d)')
        plt.ylabel('Carbohydrate intake (KJ/d)')
    else: 
        # Axes labelled in illustrator (or similar) as a panel in a larger plot
        plt.xlabel('')
        plt.ylabel('')
    # Ensure graphs always plotted on the same axis ranges. 
    # Because the markers can differ in size (transparent vs highlighted), the graph boundaries can shift. 
    plt.xlim([0, 28])
    plt.ylim([0, 38])

    plt.savefig('cage-locations_grey.png', dpi=300, transparent=True)
    plt.close()

# Individual cages
if True:
    os.makedirs('individual_cages', exist_ok=True)
    for cage_id in meta_ssb_df['Mouse No']:
        print(cage_id)
        # Plot all cage locations as highly transparent markers.
        ax = sns.scatterplot(x='prot-intake-KJ/d', y='carb-intake-KJ/d', data=meta_ssb_df, 
                        hue='%P/%C/%F', palette=palette, s=150, marker='o', alpha=0.15,
                        legend=False)
        
        # Plot the select cage as an opaque marker, slightly larger. 
        cage_index = meta_ssb_df['Mouse No'] == cage_id
        cage_df = meta_ssb_df.loc[cage_index, :]
        # By specifying both hue and palette here we ensure the same colours are maintained. 
        sns.scatterplot(x='prot-intake-KJ/d', y='carb-intake-KJ/d', data=cage_df, 
                        hue='%P/%C/%F', 
                        palette=palette_black, 
                        s=350, marker='o',
                        legend=False)
        # Plot the convex hull of the data.
        for simplex in hull.simplices:
            plt.plot(meta_nd[simplex, 0], meta_nd[simplex, 1], 'k-')

        if draw_labels:
            plt.xlabel('Protein intake (KJ/d)')
            plt.ylabel('Carbohydrate intake (KJ/d)')
        else: 
            # Axes labelled in illustrator (or similar) as a panel in a larger plot
            plt.xlabel('')
            plt.ylabel('')

        # Ensure graphs always plotted on the same axis ranges. 
        # Because the markers can differ in size (transparent vs highlighted), the graph boundaries can shift. 
        plt.xlim([0, 28])
        plt.ylim([0, 38])

        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        plt.savefig('individual_cages_grey/cage-{:d}.png'.format(cage_id), dpi=300, transparent=True)
        plt.close()

# The white standard chow "average" mouse control. 
if False:
    # Plot all cage locations as highly transparent markers.
    ax = sns.scatterplot(x='prot-intake-KJ/d', y='carb-intake-KJ/d', data=meta_ssb_df, 
                    hue='%P/%C/%F', palette=palette, s=150, marker='o', alpha=0.15,
                    legend=False)
    white_indices = meta_df['Experiment'] == 'SSB-sim-control'
    white_df = meta_df.loc[white_indices, :] 
    print(white_df)
    # By specifying both hue and palette here we ensure the same colours are maintained. 
    palette['19/63/18'] = 'k'  # For the white standard chow simulated control. 
    sns.scatterplot(x='prot-intake-KJ/d', y='carb-intake-KJ/d', data=white_df, 
                    hue='%P/%C/%F', palette=palette, s=350, marker='o',
                    legend=False)
    # Plot the convex hull of the data.
    for simplex in hull.simplices:
        plt.plot(meta_nd[simplex, 0], meta_nd[simplex, 1], 'k-')

    if draw_labels:
        plt.xlabel('Protein intake (KJ/d)')
        plt.ylabel('Carbohydrate intake (KJ/d)')
    else: 
        # Axes labelled in illustrator (or similar) as a panel in a larger plot
        plt.xlabel('')
        plt.ylabel('')

    # Ensure graphs always plotted on the same axis ranges. 
    # Because the markers can differ in size (transparent vs highlighted), the graph boundaries can shift. 
    plt.xlim([0, 28])
    plt.ylim([0, 38])

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.savefig('individual_cages_grey/cage-white_chow.png', dpi=300, transparent=True)
    plt.close()
print(meta_ssb_df.columns)

