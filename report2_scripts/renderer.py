from __future__ import annotations

import ctypes
import os

import glfw
import numpy as np
from OpenGL.GL import *
from OpenGL.GL import shaders


def _look_at(eye: np.ndarray, center: np.ndarray, up: np.ndarray) -> np.ndarray:
    """Compute a standard OpenGL look-at view matrix."""
    f = center - eye
    f = f / np.linalg.norm(f)
    s = np.cross(f, up)
    s = s / np.linalg.norm(s)
    u = np.cross(s, f)
    M = np.eye(4, dtype=np.float32)
    M[0, :3] = s
    M[1, :3] = u
    M[2, :3] = -f
    M[0, 3] = -np.dot(s, eye)
    M[1, 3] = -np.dot(u, eye)
    M[2, 3] = np.dot(f, eye)
    return M


def _perspective(fov: float, aspect: float, near: float, far: float) -> np.ndarray:
    """Compute a symmetric perspective projection matrix."""
    f = 1.0 / np.tan(np.radians(fov) / 2.0)
    M = np.zeros((4, 4), dtype=np.float32)
    M[0, 0] = f / aspect
    M[1, 1] = f
    M[2, 2] = (far + near) / (near - far)
    M[2, 3] = (2.0 * far * near) / (near - far)
    M[3, 2] = -1.0
    return M


def init_opengl_context(width: int = 1400, height: int = 900,
                        title: str = "VQE Molecular Explorer") -> object:
    """Initialise a GLFW window with an OpenGL 3.3 core-profile context.

    Steps:
      1. glfw.init() — raise RuntimeError on failure.
      2. Set window hints: CONTEXT_VERSION_MAJOR=3, MINOR=3,
         OPENGL_PROFILE=CORE_PROFILE, SAMPLES=4 (MSAA).
      3. Create window, make context current.
      4. Enable: GL_DEPTH_TEST, GL_MULTISAMPLE, GL_BLEND
         with glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA).

    Args:
        width: Window width in pixels.
        height: Window height in pixels.
        title: Window title string.

    Returns:
        The GLFW window handle.
    """
    if not glfw.init():
        raise RuntimeError("Failed to initialise GLFW")
    
    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    glfw.window_hint(glfw.SAMPLES, 4)  # Enable 4x MSAA
    
    window = glfw.create_window(width, height, title, None, None)
    if not window:
        glfw.terminate()
        raise RuntimeError("Failed to create GLFW window")
    
    glfw.make_context_current(window)
    glfw.swap_interval(1)  # VSync
    
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_MULTISAMPLE)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    
    return window


def compile_shader_program(vertex_source: str, fragment_source: str) -> int:
    """Compile and link a vertex + fragment shader into a program."""
    vertex_shader = shaders.compileShader(vertex_source, GL_VERTEX_SHADER)
    fragment_shader = shaders.compileShader(fragment_source, GL_FRAGMENT_SHADER)
    shader_program = shaders.compileProgram(vertex_shader, fragment_shader)
    return shader_program


def load_shader_from_file(shader_dir: str, vert_name: str, frag_name: str) -> int:
    """Load and compile shaders from files."""
    vert_path = os.path.join(shader_dir, vert_name)
    frag_path = os.path.join(shader_dir, frag_name)
    
    with open(vert_path, 'r') as f:
        vert_source = f.read()
    with open(frag_path, 'r') as f:
        frag_source = f.read()
    
    return compile_shader_program(vert_source, frag_source)

def create_orbital_vao(mesh: object) -> tuple[int, int, int, int]:
    """Create an OpenGL VAO/VBO/EBO for an isosurface mesh.

    Vertex layout: [x, y, z, nx, ny, nz, sign] — 7 floats, stride = 28 bytes.

    Pack mesh.vertices, mesh.normals, and mesh.signs into a single
    (n_verts, 7) float32 array. Upload to VBO via glBufferData (GL_DYNAMIC_DRAW).
    Upload mesh.faces to EBO via glBufferData (GL_DYNAMIC_DRAW).

    Vertex attribute pointers:
      - Location 0: position — 3 floats, offset 0
      - Location 1: normal — 3 floats, offset 12
      - Location 2: sign — 1 float, offset 24

    Args:
        mesh: IsosurfaceMesh with vertices, normals, faces, and signs.

    Returns:
        Tuple of (vao, vbo, ebo, n_indices) where n_indices = len(faces) * 3.
    """
    # Pack vertex data: [x, y, z, nx, ny, nz, sign]
    vertex_data = np.hstack([
        mesh.vertices, 
        mesh.normals, 
        mesh.signs[:, np.newaxis]
    ]).astype(np.float32)
    
    # Create VAO, VBO, EBO
    vao = glGenVertexArrays(1)
    vbo = glGenBuffers(1)
    ebo = glGenBuffers(1)
    
    glBindVertexArray(vao)
    
    # Upload vertex data
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertex_data.nbytes, vertex_data.flatten(), GL_DYNAMIC_DRAW)
    
    # Upload index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.faces.nbytes, mesh.faces.flatten(), GL_DYNAMIC_DRAW)
    
    # Set up vertex attributes
    stride = 7 * 4  # 7 floats * 4 bytes
    
    # Position (location 0)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    
    # Normal (location 1)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)
    
    # Sign (location 2)
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(24))
    glEnableVertexAttribArray(2)
    
    glBindVertexArray(0)
    
    n_indices = len(mesh.faces) * 3
    return vao, vbo, ebo, n_indices
    


def update_orbital_vao(vbo: int, ebo: int, mesh: object) -> int:
    """Update existing VBO/EBO with new mesh data (buffer orphaning).

    Re-pack vertex data in same [x,y,z,nx,ny,nz,sign] layout and
    re-upload via glBufferData (GL_DYNAMIC_DRAW), which orphans the
    old buffer for efficient double-buffered updates.

    Args:
        vbo: Existing vertex buffer object ID.
        ebo: Existing element buffer object ID.
        mesh: New IsosurfaceMesh data.

    Returns:
        Updated n_indices (len(faces) * 3).
    """
    vertex_data = np.hstack([mesh.vertices, mesh.normals, mesh.signs[:, np.newaxis]]).astype(np.float32)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertex_data.nbytes, vertex_data.flatten(), GL_DYNAMIC_DRAW)

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.faces.nbytes, mesh.faces.flatten(), GL_DYNAMIC_DRAW)

    return len(mesh.faces) * 3


# Global unit sphere/cylinder VAO (created once)
_unit_sphere_vao = None
_unit_sphere_n_indices = 0
_unit_cylinder_vao = None
_unit_cylinder_n_indices = 0


def _create_unit_sphere(segments: int = 16, rings: int = 16) -> tuple[int, int]:
    """Create a unit sphere VAO centered at origin with radius 1."""
    global _unit_sphere_vao, _unit_sphere_n_indices
    
    vertices = []
    indices = []
    
    for i in range(rings + 1):
        phi = np.pi * i / rings
        for j in range(segments + 1):
            theta = 2.0 * np.pi * j / segments
            x = np.sin(phi) * np.cos(theta)
            y = np.cos(phi)
            z = np.sin(phi) * np.sin(theta)
            # Position and normal are the same for a unit sphere
            vertices.extend([x, y, z, x, y, z])
    
    for i in range(rings):
        for j in range(segments):
            a = i * (segments + 1) + j
            b = a + segments + 1
            indices.extend([a, b, a + 1, b, b + 1, a + 1])
    
    vertices = np.array(vertices, dtype=np.float32)
    indices = np.array(indices, dtype=np.uint32)
    
    vao = glGenVertexArrays(1)
    vbo = glGenBuffers(1)
    ebo = glGenBuffers(1)
    
    glBindVertexArray(vao)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_STATIC_DRAW)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_STATIC_DRAW)
    
    stride = 6 * 4  # 6 floats * 4 bytes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)
    
    glBindVertexArray(0)
    
    _unit_sphere_vao = vao
    _unit_sphere_n_indices = len(indices)
    return vao, len(indices)


def _create_unit_cylinder(segments: int = 16) -> tuple[int, int]:
    """Create a unit cylinder VAO from (0,0,0) to (0,1,0) with radius 1."""
    global _unit_cylinder_vao, _unit_cylinder_n_indices
    
    vertices = []
    indices = []
    
    # Generate vertices for cylinder body
    for i in range(2):  # bottom (0) and top (1)
        y = float(i)
        for j in range(segments + 1):
            theta = 2.0 * np.pi * j / segments
            x = np.cos(theta)
            z = np.sin(theta)
            # Position and normal (normal points outward)
            vertices.extend([x, y, z, x, 0, z])
    
    # Generate indices for cylinder body
    for j in range(segments):
        a = j
        b = j + segments + 1
        indices.extend([a, b, a + 1, b, b + 1, a + 1])
    
    vertices = np.array(vertices, dtype=np.float32)
    indices = np.array(indices, dtype=np.uint32)
    
    vao = glGenVertexArrays(1)
    vbo = glGenBuffers(1)
    ebo = glGenBuffers(1)
    
    glBindVertexArray(vao)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_STATIC_DRAW)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_STATIC_DRAW)
    
    stride = 6 * 4
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)
    
    glBindVertexArray(0)
    
    _unit_cylinder_vao = vao
    _unit_cylinder_n_indices = len(indices)
    return vao, len(indices)


def draw_unit_sphere() -> None:
    """Draw the unit sphere using its VAO."""
    global _unit_sphere_vao, _unit_sphere_n_indices
    if _unit_sphere_vao is None:
        _create_unit_sphere()
    glBindVertexArray(_unit_sphere_vao)
    glDrawElements(GL_TRIANGLES, _unit_sphere_n_indices, GL_UNSIGNED_INT, None)
    glBindVertexArray(0)


def draw_unit_cylinder() -> None:
    """Draw the unit cylinder using its VAO."""
    global _unit_cylinder_vao, _unit_cylinder_n_indices
    if _unit_cylinder_vao is None:
        _create_unit_cylinder()
    glBindVertexArray(_unit_cylinder_vao)
    glDrawElements(GL_TRIANGLES, _unit_cylinder_n_indices, GL_UNSIGNED_INT, None)
    glBindVertexArray(0)


def build_atom_model_matrix(position: np.ndarray, radius: float) -> np.ndarray:
    """Build a 4x4 model matrix for an atom sphere (translate + scale)."""
    M = np.eye(4, dtype=np.float32)
    M[0, 0] = radius
    M[1, 1] = radius
    M[2, 2] = radius
    M[0, 3] = position[0]
    M[1, 3] = position[1]
    M[2, 3] = position[2]
    return M


def build_bond_model_matrix(pos_i: np.ndarray, pos_j: np.ndarray, 
                            radius: float = 0.1) -> np.ndarray:
    """Build a 4x4 model matrix for a bond cylinder between two atoms."""
    direction = pos_j - pos_i
    length = np.linalg.norm(direction)
    if length < 1e-6:
        return np.eye(4, dtype=np.float32)
    
    direction = direction / length
    
    # Default cylinder is along Y-axis, we need to rotate to direction
    up = np.array([0.0, 1.0, 0.0])
    
    # Handle case where direction is parallel to up
    if abs(np.dot(direction, up)) > 0.999:
        right = np.array([1.0, 0.0, 0.0])
    else:
        right = np.cross(up, direction)
        right = right / np.linalg.norm(right)
    
    new_up = np.cross(direction, right)
    
    # Build rotation matrix
    R = np.eye(4, dtype=np.float32)
    R[0, :3] = right
    R[1, :3] = direction
    R[2, :3] = new_up
    R = R.T  # Transpose for column-major
    
    # Scale matrix
    S = np.eye(4, dtype=np.float32)
    S[0, 0] = radius
    S[1, 1] = length
    S[2, 2] = radius
    
    # Translation matrix
    T = np.eye(4, dtype=np.float32)
    T[0, 3] = pos_i[0]
    T[1, 3] = pos_i[1]
    T[2, 3] = pos_i[2]
    
    return T @ R @ S


def draw_molecule(atom_positions: np.ndarray, atom_symbols: list[str],
                  bonds: list[tuple[int, int]], shader_program: int) -> None:
    """Render atoms as spheres and bonds as cylinders using instanced rendering.

    Atom colours by element:
      H=white(1,1,1), C=grey(0.5,0.5,0.5), N=blue(0,0,1),
      O=red(1,0,0), Li=purple(0.6,0.2,0.8), F=green(0,1,0).

    Van der Waals radii (scaled by 0.3): H=1.2, C=1.7, N=1.55, O=1.52, Li=1.82.

    For each atom: build model matrix (translate + scale by radius),
    set shader uniforms, draw unit sphere VAO.

    For each bond (i,j): build model matrix transforming a unit cylinder
    from (0,0,0)->(0,0,1) to atom_positions[i]->atom_positions[j]
    with radius=0.1, set shader uniforms, draw unit cylinder VAO.

    Args:
        atom_positions: Array of atom positions, shape (n_atoms, 3).
        atom_symbols: List of element symbols.
        bonds: List of (atom_i, atom_j) index pairs.
        shader_program: Compiled OpenGL shader program ID.
    """
    element_colors = {
        "H": (1.0, 1.0, 1.0),
        "C": (0.5, 0.5, 0.5),
        "N": (0.0, 0.0, 1.0),
        "O": (1.0, 0.0, 0.0),
        "Li": (0.6, 0.2, 0.8),
        "F": (0.0, 1.0, 0.0)
    }
    vdw_radii = {
        "H": 1.2,
        "C": 1.7,
        "N": 1.55,
        "O": 1.52,
        "Li": 1.82,
        "F": 1.47
    }
    for i, (pos, sym) in enumerate(zip(atom_positions, atom_symbols)):
        color = element_colors.get(sym, (0.5, 0.5, 0.5))  # Default grey
        radius = vdw_radii.get(sym, 1.5) * 0.3  # Default radius
        model_matrix = build_atom_model_matrix(pos, radius)
        glUniformMatrix4fv(glGetUniformLocation(shader_program, "model"), 1, GL_FALSE, model_matrix.flatten())
        glUniform3f(glGetUniformLocation(shader_program, "color"), *color)
        draw_unit_sphere()
    for i, j in bonds:
        pos_i, pos_j = atom_positions[i], atom_positions[j]
        model_matrix = build_bond_model_matrix(pos_i, pos_j, radius=0.1)
        glUniformMatrix4fv(glGetUniformLocation(shader_program, "model"), 1, GL_FALSE, model_matrix.flatten())
        glUniform3f(glGetUniformLocation(shader_program, "color"), 0.8, 0.8, 0.8)  # Bond color
        draw_unit_cylinder()




class ArcballCamera:
    """Orbiting camera with mouse-drag rotation and scroll zoom.

    Attributes:
        position: Camera position in world space, shape (3,).
        target: Look-at target point, shape (3,).
        up: Up vector, default [0, 1, 0].
        yaw: Horizontal angle in degrees (default -90).
        pitch: Vertical angle in degrees (default 20, clamped ±89).
        zoom: Distance from target (initialised from |position - target|).
    """
    position: np.ndarray
    target: np.ndarray
    up: np.ndarray
    yaw: float
    pitch: float
    zoom: float

    def __init__(self, position: np.ndarray | None = None,
                 target: np.ndarray | None = None) -> None:
        """Initialise camera.

        Defaults: position=[0,0,10], target=[0,0,0], up=[0,1,0],
        yaw=-90, pitch=20, zoom = |position - target|.

        Args:
            position: Initial camera position (or None for default).
            target: Initial look-at target (or None for default).
        """
        self.position = position if position is not None else np.array([0.0, 0.0, 10.0])
        self.target = target if target is not None else np.array([0.0, 0.0, 0.0])
        self.up = np.array([0.0, 1.0, 0.0])
        self.yaw = -90.0
        self.pitch = 20.0
        self.zoom = np.linalg.norm(self.position - self.target)

    def on_mouse_drag(self, dx: float, dy: float) -> None:
        """Update yaw/pitch from mouse drag deltas.

        yaw += dx * 0.3; pitch -= dy * 0.3 (clamped to [-89, 89]).
        Calls _update() to recompute position.

        Args:
            dx: Horizontal mouse delta in pixels.
            dy: Vertical mouse delta in pixels.
        """
        self.yaw += dx * 0.3
        self.pitch -= dy * 0.3
        self.pitch = np.clip(self.pitch, -89, 89)
        self._update()

    def on_scroll(self, delta: float) -> None:
        """Update zoom from scroll wheel delta.

        zoom -= delta * 0.5 (clamped to min 1.0).
        Calls _update() to recompute position.

        Args:
            delta: Scroll wheel delta.
        """
        self.zoom -= delta * 0.5
        self.zoom = np.clip(self.zoom, 1.0, np.inf)
        self._update()

    def _update(self) -> None:
        """Recompute camera position from yaw, pitch, zoom, and target.

        position = target + zoom * [cos(pitch)*cos(yaw), sin(pitch), cos(pitch)*sin(yaw)]
        (yaw and pitch converted to radians).
        """
        yaw_rad = np.radians(self.yaw)
        pitch_rad = np.radians(self.pitch)
        direction = np.array([
            np.cos(pitch_rad) * np.cos(yaw_rad),
            np.sin(pitch_rad),
            np.cos(pitch_rad) * np.sin(yaw_rad)
        ])
        self.position = self.target + self.zoom * direction

    def get_view_matrix(self) -> np.ndarray:
        """Compute a 4x4 look-at view matrix.

        Uses position, target, and up vector via a _look_at helper.

        Returns:
            4x4 numpy array (view matrix).
        """
        return _look_at(self.position, self.target, self.up)

    def get_projection_matrix(self, aspect: float,
                              fov: float = 45.0) -> np.ndarray:
        """Compute a 4x4 perspective projection matrix.

        Uses a _perspective helper with fov, aspect, near=0.1, far=200.0.

        Args:
            aspect: Viewport width / height.
            fov: Vertical field of view in degrees (default 45).

        Returns:
            4x4 numpy array (projection matrix).
        """
        return _perspective(fov, aspect, 0.1, 200.0)
    


