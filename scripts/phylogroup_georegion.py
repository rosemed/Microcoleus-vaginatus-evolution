import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from scipy.spatial import ConvexHull

# 设置中文字体
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 读取地理位置数据
def read_location_data(file_path):
    """读取包含菌株经纬度信息的文件"""
    try:
        # 假设文件是制表符分隔的，包含columns: strain_id, longitude, latitude, lineage
        df = pd.read_table(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return None

# 计算每个系群的地理范围
def calculate_geographic_range(df):
    """计算每个系群的中心点、边界和地理范围"""
    lineages = df['lineage'].unique()
    lineage_data = {}
    
    for lineage in lineages:
        # 过滤当前系群的数据
        lineage_df = df[df['lineage'] == lineage]
        
        # 计算中心点 (平均经纬度)
        center_lon = lineage_df['longitude'].mean()
        center_lat = lineage_df['latitude'].mean()
        
        # 计算边界 (使用凸包算法)
        points = np.array(lineage_df[['longitude', 'latitude']])
        if len(points) >= 3:  # 至少需要3个点来形成凸包
            hull = ConvexHull(points)
            boundary = points[hull.vertices]
        else:
            boundary = points  # 如果点太少，直接使用这些点作为边界
            
        # 计算地理范围 (使用边界点的经纬度范围)
        min_lon, max_lon = points[:, 0].min(), points[:, 0].max()
        min_lat, max_lat = points[:, 1].min(), points[:, 1].max()
        area = (max_lon - min_lon) * (max_lat - min_lat)
        
        lineage_data[lineage] = {
            'center': (center_lon, center_lat),
            'boundary': boundary,
            'area': area,
            'points': points
        }
    
    return lineage_data

# 可视化不同系群的地理分布和范围
def visualize_geographic_ranges(lineage_data, output_file=None):
    """可视化不同系群的地理分布和范围"""
    plt.figure(figsize=(30, 10))
    
    # 创建基础地图
    m = Basemap(
        projection='merc',
        llcrnrlon=-110, llcrnrlat=-25,
        urcrnrlon=120, urcrnrlat=60,
        resolution='l'
    )
    
    # 绘制海岸线、国家和州界
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.3)
    
    # 设置颜色映射
    colors = cm.rainbow(np.linspace(0, 1, len(lineage_data)))
    cmap = ListedColormap(colors)
    
    # 为每个系群绘制地理范围和中心点
    for i, (lineage, data) in enumerate(lineage_data.items()):
        color = colors[i]
        
        # 绘制系群的所有点
        x, y = m(data['points'][:, 0], data['points'][:, 1])
        m.scatter(x, y, s=15, c=[color], label=f'系群 {lineage}', alpha=0.7)
        
        # 绘制中心点
        center_x, center_y = m(*data['center'])
        m.plot(center_x, center_y, 'o', color='black', markersize=6, markeredgecolor='white')
        
        # 绘制边界 (如果有足够的点)
        if len(data['boundary']) >= 3:
            boundary_x, boundary_y = m(data['boundary'][:, 0], data['boundary'][:, 1])
            boundary_xy = list(zip(boundary_x, boundary_y))
            poly = Polygon(boundary_xy, facecolor=color, alpha=0.3, edgecolor='black', linewidth=1)
            plt.gca().add_patch(poly)
    
    # 添加图例和标题
    plt.legend(loc='upper right', framealpha=0.8)
    
    if output_file:
        plt.savefig(f"{output_file}.pdf", format='pdf', bbox_inches='tight')
        print(f"地图已保存到 {output_file}.pdf")
    
    plt.show()

# 主函数
def main():
    file_path = 'phylogroup_location.txt'  # 替换为实际文件路径
    df = read_location_data(file_path)
    
    if df is not None:
        # 筛选G2-G8系群
        target_lineages = ['G2', 'G4', 'G5', 'G7', 'G8']
        filtered_df = df[df['lineage'].isin(target_lineages)]
        
        if not filtered_df.empty:
            lineage_data = calculate_geographic_range(filtered_df)
            visualize_geographic_ranges(lineage_data, output_file='lineage_geographic_range')
        else:
            print("未找到系群的数据")

if __name__ == "__main__":
    main()
