import plotly.graph_objects as go
import numpy as np
import os

def plot_interactive_matrix(output_path, beta_matrix, pvalues_matrix, labels):
    """
    Plot an interactive connectivity matrix using Plotly.
    """
    print("\n📊 Plotting interactive connectivity matrix ...\n")

    # Prepare pvalues matrix
    pvalues_corrected = np.array(pvalues_matrix)
    significance_mask = np.zeros_like(pvalues_corrected, dtype=bool)
    upper_triangle_indices = np.triu_indices_from(pvalues_corrected, k=1)
    significance_mask[upper_triangle_indices] = pvalues_corrected[upper_triangle_indices] < 0.05

    # Prepare beta matrix
    beta_matrix = np.array(beta_matrix)
    abs_max = np.nanmax(np.abs(beta_matrix))
    zmin, zmax = -abs_max, abs_max

    # Create a full matrix of NaNs with same shape
    beta_matrix_upper = np.full_like(beta_matrix, np.nan)

    # Fill in only the upper triangle (excluding diagonal)
    i_upper, j_upper = np.triu_indices_from(beta_matrix, k=1)
    beta_matrix_upper[i_upper, j_upper] = beta_matrix[i_upper, j_upper]

    custom_colorscale = [
        [0.0, "blue"],   
        [0.5, "white"], 
        [1.0, "red"]     
    ]

    # Prepare beta heatmap
    heatmap = go.Heatmap(
        z=beta_matrix_upper,
        x=labels,
        y=labels,
        colorscale=custom_colorscale,
        zmin=zmin,
        zmax=zmax,
        colorbar=dict(title='Beta values'),
        hovertemplate="%{y} - %{x}<br>Beta value: %{z}<extra></extra>"
    )

    # Semi-transparent mask to indicate significance
    mask = np.where(significance_mask, 1.0, np.nan)  # show only where significant
    sig_overlay = go.Heatmap(
        z=mask,
        x=labels,
        y=labels,
        showscale=False,
        colorscale=[[0, 'rgba(255,255,0,0.6)'], [1, 'rgba(255,255,0,0.6)']],  # with yellow overlay
        hoverinfo='skip'
    )

    # Combine
    fig = go.Figure(data=[heatmap, sig_overlay])

    fig.update_layout(
        title=f'Interactive beta value heatmap with FDR correction',
        xaxis=dict(title='Labels (Colonnes)', tickangle=45),
        yaxis=dict(title='Labels (Lignes)', autorange='reversed'),
        width=1200,  
        height=800   
    )

    # Save in html file
    filename = f'interactive_heatmap'
    output_path = os.path.join(output_path, "{}.html".format(filename))
    fig.write_html(
    output_path,
    include_plotlyjs='cdn', 
    full_html=True,         
    config={"responsive": True} 
    )
    print(f"✅ Heatmap interactive saved in : {output_path}")