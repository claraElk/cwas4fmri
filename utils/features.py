import json

def extract_feature_settings(json_file_path):
    """
    Extract feature names and their settings from JSON file and create a table
    
    Args:
        json_file_path (str): Path to the JSON file containing feature settings
    
    Returns:
        tuple: (pd.DataFrame, str) containing feature settings and the common atlas name
    
    Raises:
        ValueError: If different atlases are found across features
    """
    print("Identify variables based on spec.json ...")
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    
    # Extract features and settings
    feature_settings = []
    feature_atlases = []
    
    for feature in data['features']:
        if feature['type'] == 'atlas_based_connectivity':
            feature_settings.append(feature['name'])
            feature_atlases.append(feature['atlases'][0])
            
    # Check if all features use the same atlas
    if len(set(feature_atlases)) > 1:
        raise ValueError(f"Error: Multiple different atlases found across features:\n{set(feature_atlases)}")
    
    # Get the common atlas name
    common_atlas = feature_atlases[0] if feature_atlases else None
    
    # Print results
    print("\nFeature Settings Found in spec.json:")
    print(feature_settings)
    
    print("\nAtlas Information:")
    if common_atlas:
        print(f"Common atlas used across all features: {common_atlas}")
    else:
        print("No atlas-based connectivity features found")
    
    return feature_settings, common_atlas