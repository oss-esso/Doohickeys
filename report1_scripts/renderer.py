"""OpenGL 3.3 real-time renderer for the QAOA energy landscape."""

from __future__ import annotations

import ctypes

import glfw
from OpenGL.GL import *
from OpenGL.GL import shaders
import numpy as np


def init_opengl_context(
    width: int = 1280,
    height: int = 720,
    title: str = "QAOA Energy Landscape Explorer",
) -> glfw._GLFWwindow:
    """Create a GLFW window with an OpenGL 3.3 core-profile context."""
    if not glfw.init():
        raise RuntimeError("Failed to initialise GLFW")
    glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    glfw.window_hint(glfw.SAMPLES, 4)
    window = glfw.create_window(width, height, title, None, None)
    if not window:
        glfw.terminate()
        raise RuntimeError("Failed to create GLFW window")
    glfw.make_context_current(window)
    glfw.swap_interval(1)
    glEnable(GL_MULTISAMPLE)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
    return window


def compile_shader_program(
    vertex_source: str,
    fragment_source: str,
) -> int:
    """Compile and link a vertex + fragment shader into a program."""
    vertex_shader = shaders.compileShader(vertex_source, GL_VERTEX_SHADER)
    fragment_shader = shaders.compileShader(fragment_source, GL_FRAGMENT_SHADER)
    shader_program = shaders.compileProgram(vertex_shader, fragment_shader)
    return shader_program


def create_surface_vbo(
    gamma_vals: np.ndarray,
    beta_vals: np.ndarray,
    energies: np.ndarray,
) -> tuple[int, int, int]:
    """Build a triangle-mesh VBO from the (γ, β, energy) grid.

    Vertex layout (7 floats per vertex, stride 28 bytes):
      [x=γ, y=β, z=energy, nx, ny, nz, energy_normalised]
    """
    nrows, ncols = energies.shape
    e_min, e_max = energies.min(), energies.max()
    e_range = e_max - e_min if e_max != e_min else 1.0

    # Compute normals via finite differences
    dEdg = np.gradient(energies, gamma_vals, axis=0)
    dEdb = np.gradient(energies, beta_vals, axis=1)

    vertices = []
    for i in range(nrows):
        for j in range(ncols):
            x, y, z = gamma_vals[i], beta_vals[j], energies[i, j]
            n = np.array([-dEdg[i, j], -dEdb[i, j], 1.0])
            n /= np.linalg.norm(n)
            e_norm = (energies[i, j] - e_min) / e_range
            vertices.extend([x, y, z, n[0], n[1], n[2], e_norm])

    vertices = np.array(vertices, dtype=np.float32)

    # Index buffer: two triangles per grid cell
    indices = []
    for i in range(nrows - 1):
        for j in range(ncols - 1):
            tl = i * ncols + j
            tr = tl + 1
            bl = (i + 1) * ncols + j
            br = bl + 1
            indices.extend([tl, bl, tr, tr, bl, br])
    indices = np.array(indices, dtype=np.uint32)

    vao = glGenVertexArrays(1)
    glBindVertexArray(vao)

    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_DYNAMIC_DRAW)

    ebo = glGenBuffers(1)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_STATIC_DRAW)

    stride = 7 * 4  # 7 floats × 4 bytes
    # location 0: position (vec3)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    # location 1: normal (vec3)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)
    # location 2: energy normalised (float)
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(24))
    glEnableVertexAttribArray(2)

    glBindVertexArray(0)
    return vao, vbo, ebo


def update_surface_vbo(
    vbo_id: int,
    gamma_vals: np.ndarray,
    beta_vals: np.ndarray,
    energies: np.ndarray,
) -> None:
    """Stream new landscape data into an existing VBO via glBufferSubData."""
    nrows, ncols = energies.shape
    e_min, e_max = energies.min(), energies.max()
    e_range = e_max - e_min if e_max != e_min else 1.0

    dEdg = np.gradient(energies, gamma_vals, axis=0)
    dEdb = np.gradient(energies, beta_vals, axis=1)

    vertices = []
    for i in range(nrows):
        for j in range(ncols):
            x, y, z = gamma_vals[i], beta_vals[j], energies[i, j]
            n = np.array([-dEdg[i, j], -dEdb[i, j], 1.0])
            n /= np.linalg.norm(n)
            e_norm = (energies[i, j] - e_min) / e_range
            vertices.extend([x, y, z, n[0], n[1], n[2], e_norm])

    vertices = np.array(vertices, dtype=np.float32)
    glBindBuffer(GL_ARRAY_BUFFER, vbo_id)
    glBufferSubData(GL_ARRAY_BUFFER, 0, vertices.nbytes, vertices)


def draw_optimizer_trajectory(
    trajectory: list[tuple[float, float, float]],
    shader_program: int,
) -> None:
    """Render the optimizer path as a GL_LINE_STRIP over the surface."""
    if len(trajectory) < 2:
        return

    n = len(trajectory)
    data = []
    for idx, (g, b, e) in enumerate(trajectory):
        t = idx / max(n - 1, 1)
        data.extend([g, b, e + 0.01, t])
    data = np.array(data, dtype=np.float32)

    vao = glGenVertexArrays(1)
    glBindVertexArray(vao)
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, data.nbytes, data, GL_STATIC_DRAW)

    stride = 4 * 4
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
    glEnableVertexAttribArray(0)
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(12))
    glEnableVertexAttribArray(1)

    glUseProgram(shader_program)
    glLineWidth(3.0)
    glDrawArrays(GL_LINE_STRIP, 0, n)

    glBindVertexArray(0)
    glDeleteBuffers(1, [vbo])
    glDeleteVertexArrays(1, [vao])


class ArcballCamera:
    """Orbit camera controlled by mouse drag (yaw/pitch) and scroll (zoom)."""

    def __init__(
        self,
        position: np.ndarray = np.array([0.0, 0.0, 5.0]),
        target: np.ndarray = np.array([0.0, 0.0, 0.0]),
    ) -> None:
        self.target = target.astype(np.float32).copy()
        self.up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
        diff = position - target
        self.zoom = float(np.linalg.norm(diff))
        self.yaw = float(np.degrees(np.arctan2(diff[2], diff[0])))
        self.pitch = float(np.degrees(np.arcsin(np.clip(diff[1] / self.zoom, -1, 1))))
        self.position = position.astype(np.float32).copy()

    def on_mouse_drag(self, dx: float, dy: float) -> None:
        sensitivity = 0.3
        self.yaw += dx * sensitivity
        self.pitch += dy * sensitivity
        self.pitch = np.clip(self.pitch, -89.0, 89.0)
        self._update_position()

    def on_scroll(self, delta: float) -> None:
        self.zoom -= delta * 0.5
        self.zoom = max(self.zoom, 0.5)
        self._update_position()

    def _update_position(self) -> None:
        pitch_rad = np.radians(self.pitch)
        yaw_rad = np.radians(self.yaw)
        self.position = self.target + self.zoom * np.array([
            np.cos(pitch_rad) * np.cos(yaw_rad),
            np.sin(pitch_rad),
            np.cos(pitch_rad) * np.sin(yaw_rad),
        ], dtype=np.float32)

    def get_view_matrix(self) -> np.ndarray:
        return _look_at(self.position, self.target, self.up)

    def get_projection_matrix(
        self,
        aspect_ratio: float,
        fov: float = 45.0,
    ) -> np.ndarray:
        return _perspective(fov, aspect_ratio, 0.1, 100.0)


def _look_at(
    eye: np.ndarray,
    center: np.ndarray,
    up: np.ndarray,
) -> np.ndarray:
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


def _perspective(
    fov: float,
    aspect: float,
    near: float,
    far: float,
) -> np.ndarray:
    """Compute a symmetric perspective projection matrix."""
    f = 1.0 / np.tan(np.radians(fov) / 2.0)
    M = np.zeros((4, 4), dtype=np.float32)
    M[0, 0] = f / aspect
    M[1, 1] = f
    M[2, 2] = (far + near) / (near - far)
    M[2, 3] = (2.0 * far * near) / (near - far)
    M[3, 2] = -1.0
    return M
