#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in float aProgress;  // 0..1

uniform mat4 uMVP;

out float vProgress;

void main() {
    vProgress = aProgress;
    gl_Position = uMVP * vec4(aPosition, 1.0);
}
