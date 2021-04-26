#version 330

uniform float radius;
uniform mat4 projection_matrix;

in block
{
	flat vec4 mv_pos;
	flat uint level;
	flat uint isBorder;
}
In;

out vec4 out_color;

const uint NUM_COLORS = 4u;

const vec3 levelColors[NUM_COLORS] =
vec3[NUM_COLORS](
    vec3(1.0, 0.0, 0.0),
    vec3(1.0, 1.0, 0.0),
    vec3(0.0, 0.0, 1.0),
    vec3(0.3, 0.3, 0.3) // Region border color
);

void main()
{
    // calculate normal
    vec3 n;
    n.xy = gl_PointCoord* 2.0 - vec2(1.0);
    float mag = dot(n.xy, n.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
    n.z = sqrt(1.0 - mag);

	vec3 eye = In.mv_pos.xyz + vec3(0.0, 0.0, radius * n.z);
	float depth = (projection_matrix[2][2] * eye.z + projection_matrix[3][2])
        / (projection_matrix[2][3] * eye.z + projection_matrix[3][3]);

    gl_FragDepth = (depth + 1.0) / 2.0;

	// compute final color
    uint colorIndex = In.level * (1u - In.isBorder) + (NUM_COLORS - 1u) * In.isBorder;
    out_color = vec4(levelColors[colorIndex], 1.0);
}
