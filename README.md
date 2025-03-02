# TVVViewer - 3D Gadget Viewer

TVVViewer is a web-based 3D visualization tool built with [trame](https://github.com/Kitware/trame) and [VTK](https://vtk.org/). It allows you to interactively view and manipulate VTK files organized in subdirectories. The application leverages Vue2 and Vuetify for its UI and offers features like adjustable representations, color mapping, vector warping, orientation markers, and cube axes.

## Features

- **Interactive Visualization:** Load and explore VTK files with various rendering options.
- **Multiple Representations:** Switch between points, wireframe, surface, and surface with edges.
- **Color Mapping:** Choose from different color lookup tables including Rainbow, Inverted Rainbow, Greyscale, and Inverted Greyscale.
- **Vector Warping:** Apply a warp vector filter to emphasize vector data.
- **UI Controls:** Easily adjust visualization parameters via a modern web interface.
- **Annotations and Markers:** Includes a scalar bar, cube axes, and an orientation marker for better spatial understanding.

## Requirements

- **Python 3.9**
- **pip**

## Installation

1. **Upgrade pip:**

   ```bash
   python -m pip install --upgrade pip
   ```

2. **Install Dependencies:**

   ```bash
   pip install "trame>=2.4.2"
   pip install "vtk>=9.1.0"
   pip install -e .
   pip install trame-vtk
   pip install trame-vuetify
   ```

## Directory Structure

```
.
├── app.py              # Main application file
├── data/               # Directory containing subdirectories with VTK files
├── README.md           # This file
└── ...                 # Other project files and resources
```

## Usage

After installing the dependencies, run the application with:

```bash
python app.py
```

This command starts the server and launches a web interface where you can:
- Select VTK files from different subdirectories within the `./data` folder.
- Adjust visualization parameters such as warp factor, mesh opacity, and representation mode.
- Change the color mapping and view additional markers like the cube axes and orientation marker.

## Code Overview

The program is organized as follows:

- **VTK Pipeline:**  
  Reads and processes VTK files using VTK readers and filters. It applies operations like vector warping and scalar extraction, and configures mappers and actors for rendering.

- **Trame Setup:**  
  Initializes the trame server and sets up the Vue2 client. The application uses a single-page layout with a drawer panel for various UI controls.

- **UI Components:**  
  Uses Vuetify components to create interactive controls for adjusting visualization parameters (e.g., sliders, selectors, and buttons). These controls update the VTK visualization in real time.

- **State Management:**  
  Uses trame's state mechanism to handle changes in visualization settings (such as background color, scale factor, representation mode, and file selection).
