#version 330 core

in vec3 vFragPos;
in vec3 vNormal;
in float vSign;

uniform vec3 uLightPos;
uniform vec3 uViewPos;

out vec4 FragColor;

void main() {
    // Lobe colouring: positive = blue, negative = red
    vec3 posColor = vec3(0.2, 0.4, 0.9);   // blue
    vec3 negColor = vec3(0.9, 0.2, 0.2);   // red
    vec3 baseColor = mix(negColor, posColor, (vSign + 1.0) * 0.5);

    // Phong lighting
    vec3 norm = normalize(vNormal);
    vec3 lightDir = normalize(uLightPos - vFragPos);
    vec3 viewDir = normalize(uViewPos - vFragPos);
    vec3 reflectDir = reflect(-lightDir, norm);

    // Two-sided lighting (orbital lobes seen from inside)
    float diff = abs(dot(norm, lightDir));

    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * baseColor;
    vec3 diffuse = diff * baseColor;

    float specularStrength = 0.4;
    float spec = pow(max(abs(dot(viewDir, reflectDir)), 0.0), 32.0);
    vec3 specular = specularStrength * spec * vec3(1.0);

    vec3 result = ambient + diffuse + specular;
    FragColor = vec4(result, 0.7);  // semi-transparent
}
