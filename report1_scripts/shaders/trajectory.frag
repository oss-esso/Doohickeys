#version 330 core
in float vProgress;
out vec4 FragColor;

void main() {
    // Blue (start) → Red (end)
    vec3 color = mix(vec3(0.0, 0.3, 1.0), vec3(1.0, 0.1, 0.0), vProgress);
    FragColor = vec4(color, 1.0);
}
