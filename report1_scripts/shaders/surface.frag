#version 330 core

in vec3 vFragPos;
in vec3 vNormal;
in float vEnergy;

uniform vec3 uLightPos;
uniform vec3 uViewPos;

out vec4 FragColor;

// HSV to RGB conversion
vec3 hsv2rgb(vec3 c) {
    vec3 p = abs(fract(c.xxx + vec3(1.0, 2.0/3.0, 1.0/3.0)) * 6.0 - 3.0);
    return c.z * mix(vec3(1.0), clamp(p - 1.0, 0.0, 1.0), c.y);
}

void main() {
    // Energy → Hue mapping: low energy = blue (0.66), high = red (0.0)
    float hue = (1.0 - vEnergy) * 0.66;
    vec3 baseColor = hsv2rgb(vec3(hue, 0.9, 0.9));
    
    // Phong lighting
    vec3 norm = normalize(vNormal);
    vec3 lightDir = normalize(uLightPos - vFragPos);
    vec3 viewDir = normalize(uViewPos - vFragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    
    // Ambient
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * baseColor;
    
    // Diffuse
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * baseColor;
    
    // Specular
    float specularStrength = 0.5;
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = specularStrength * spec * vec3(1.0);
    
    vec3 result = ambient + diffuse + specular;
    FragColor = vec4(result, 0.9);  // slight transparency
}
