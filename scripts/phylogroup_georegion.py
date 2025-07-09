import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from scipy.spatial import ConvexHull

plt.rcParams["font.family"] = ["Arial"]
plt.rcParams['axes.unicode_minus'] = False  

# Read geographical location data
def read_location_data(file_path):
    """Read the file containing the longitude and latitude information of strains"""
    try:
        # Assume the file is tab-separated and contains columns: strain_id, longitude, latitude, phylogroup
        df = pd.read_table(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"Error occurred while reading the file: {e}")
        return None

# Calculate the geographical range of each phylogroup
def calculate_geographic_range(df):
    """Calculate the center point, boundary, and geographical range of each phylogroup"""
    phylogroups = df['phylogroup'].unique()
    phylogroup_data = {}
    
    for phylogroup in phylogroups:
        # Filter the data of the current phylogroup
        phylogroup_df = df[df['phylogroup'] == phylogroup]
        
        # Calculate the center point (average longitude and latitude)
        center_lon = phylogroup_df['longitude'].mean()
        center_lat = phylogroup_df['latitude'].mean()
        
        # Calculate the boundary (using the convex hull algorithm)
        points = np.array(phylogroup_df[['longitude', 'latitude']])
        if len(points) >= 3:  # At least 3 points are required to form a convex hull
            hull = ConvexHull(points)
            boundary = points[hull.vertices]
        else:
            boundary = points  # If there are too few points, directly use these points as the boundary
            
        # Calculate the geographical range (using the longitude and latitude range of the boundary points)
        min_lon, max_lon = points[:, 0].min(), points[:, 0].max()
        min_lat, max_lat = points[:, 1].min(), points[:, 1].max()
        area = (max_lon - min_lon) * (max_lat - min_lat)
        
        phylogroup_data[phylogroup] = {
            'center': (center_lon, center_lat),
            'boundary': boundary,
            'area': area,
            'points': points
        }
    
    return phylogroup_data

# Visualize the geographical distribution and range of different phylogroups
def visualize_geographic_ranges(phylogroup_data, output_file=None):
    """Visualize the geographical distribution and range of different phylogroups"""
    plt.figure(figsize=(30, 10))
    
    # Create a base map
    m = Basemap(
        projection='merc',
        llcrnrlon=-110, llcrnrlat=-25,
        urcrnrlon=120, urcrnrlat=60,
        resolution='l'
    )
    
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.3)
    
    # Set the color map
    colors = cm.rainbow(np.linspace(0, 1, len(phylogroup_data)))
    cmap = ListedColormap(colors)
    
    # Draw the geographical range and center point for each phylogroup
    for i, (phylogroup, data) in enumerate(phylogroup_data.items()):
        color = colors[i]
        
        # Draw all points of the phylogroup
        x, y = m(data['points'][:, 0], data['points'][:, 1])
        m.scatter(x, y, s=15, c=[color], label=f'系群 {phylogroup}', alpha=0.7)
        
        # Draw the center point
        center_x, center_y = m(*data['center'])
        m.plot(center_x, center_y, 'o', color='black', markersize=6, markeredgecolor='white')
        
        # Draw the boundary (if there are enough points)
        if len(data['boundary']) >= 3:
            boundary_x, boundary_y = m(data['boundary'][:, 0], data['boundary'][:, 1])
            boundary_xy = list(zip(boundary_x, boundary_y))
            poly = Polygon(boundary_xy, facecolor=color, alpha=0.3, edgecolor='black', linewidth=1)
            plt.gca().add_patch(poly)
    
    # Add a legend and title
    plt.legend(loc='upper right', framealpha=0.8)
    
    if output_file:
        plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight')
        print(f"The map has been saved to {output_file}.pdf")
    
    plt.show()

def main():
    file_path = 'phylogroup_location.txt'  # Replace with the actual file path
    df = read_location_data(file_path)
    
    if df is not None:
        # Filter phylogroups
        target_phylogroups = ['G2', 'G4', 'G5', 'G7', 'G8']
        filtered_df = df[df['phylogroup'].isin(target_phylogroups)]
        
        if not filtered_df.empty:
            phylogroup_data = calculate_geographic_range(filtered_df)
            visualize_geographic_ranges(phylogroup_data, output_file='phylogroup_geographic_range')
        else:
            print("No data found for the phylogroups")

if __name__ == "__main__":
    main()
