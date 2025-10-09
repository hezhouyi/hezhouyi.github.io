import os

def generate_gallery_html():
    # Path to the images folder
    image_folder = "image"
    
    # Get all image files from the folder
    image_extensions = ('.jpg', '.jpeg', '.png', '.gif', '.bmp', '.webp', '.JPG', '.JPEG', '.PNG', '.GIF', '.BMP', '.WEBP','.heic')
    
    if not os.path.exists(image_folder):
        print(f"Error: '{image_folder}' folder not found. Please create it and add images.")
        return
    
    images = [f for f in os.listdir(image_folder) if f.endswith(image_extensions)]
    
    if not images:
        print(f"No images found in '{image_folder}' folder.")
        return
    
    # Sort images alphabetically
    images.sort()
    
    # Generate HTML content
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Joey's Photo Gallery ~</title>
    <link rel="stylesheet" href="style.css">
</head>
<body>
    <h1>Welcome to Joey's Photo Gallery</h1>
    <p>Click on any image to view it in full size, hope you will enjoy ~(d＇∀＇d)~ </p>
"""
    
    # Generate gallery items for each image
    for image in images:
        # Get image name without extension for description
        image_name = os.path.splitext(image)[0]
        
        html_content += f"""    <div class="gallery">
        <a target="_blank" href="{image_folder}/{image}">
            <img src="{image_folder}/{image}" alt="{image_name}">
        </a>
        <div class="description">{image_name}</div>
    </div>

"""
    
    html_content += """</body>
</html>"""
    
    # Write to photography.html
    with open("photography.html", "w") as f:
        f.write(html_content)
    
    print(f"✓ Generated photography.html with {len(images)} images")
    print(f"Images included: {', '.join(images)}")

if __name__ == "__main__":
    generate_gallery_html()
