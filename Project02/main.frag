// ~ 100 spheres, 10 samples, using Schlick Approximation

#define MAX_BOUNCES 10
#define MIN_DIST 0.001
#define MAX_DIST 100.0

struct Material {
    vec3 albedo;
    float fuzz;
    float ref_idx;
    int type; // 0: lambertian, 1: metal, 2: dielectric
};

struct HitRecord {
    float t;
    vec3 p;
    vec3 normal;
    Material material;
    bool hit;
};

// Hash function for pseudo-random numbers
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453123);
}

vec3 random_in_unit_sphere(vec2 seed) {
    float z = hash(seed) * 2.0 - 1.0;
    float a = hash(seed + 1.3) * 6.2831853;
    float r = sqrt(1.0 - z*z);
    return vec3(r * cos(a), r * sin(a), z);
}


bool hit_sphere(vec3 center, float radius, vec3 ro, vec3 rd, out float t, out vec3 normal) {
    vec3 oc = ro - center;
    float a = dot(rd, rd);
    float b = dot(oc, rd);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0.0) {
        float sqrtd = sqrt(discriminant);
        float temp = (-b - sqrtd) / a;
        if (temp > MIN_DIST && temp < MAX_DIST) {
            t = temp;
            vec3 hit = ro + rd * t;
            normal = (hit - center) / radius;
            return true;
        }
        temp = (-b + sqrtd) / a;
        if (temp > MIN_DIST && temp < MAX_DIST) {
            t = temp;
            vec3 hit = ro + rd * t;
            normal = (hit - center) / radius;
            return true;
        }
    }
    return false;
}

HitRecord world_hit(vec3 ro, vec3 rd) {
    HitRecord rec;
    rec.t = MAX_DIST;
    rec.hit = false;

    // Ground
    float t;
    vec3 normal;
    if (hit_sphere(vec3(0, -1000, 0), 1000.0, ro, rd, t, normal)) {
        rec.t = t;
        rec.p = ro + t * rd;
        rec.normal = normal;
        rec.material = Material(vec3(0.5), 0.0, 1.0, 0);
        rec.hit = true;
    }

    // Main spheres
    vec3 centers[3] = vec3[3](vec3(0,1,0), vec3(-4,1,0), vec3(4,1,0));
    Material mats[3] = Material[3](
        Material(vec3(1.0), 0.0, 1.5, 2),
        Material(vec3(0.4, 0.2, 0.1), 0.0, 1.0, 0),
        Material(vec3(0.7, 0.6, 0.5), 0.0, 1.0, 1)
    );

    for (int i = 0; i < 3; i++) {
        if (hit_sphere(centers[i], 1.0, ro, rd, t, normal) && t < rec.t) {
            rec.t = t;
            rec.p = ro + t * rd;
            rec.normal = normal;
            rec.material = mats[i];
            rec.hit = true;
        }
    }

    // Procedural small spheres
    for (int a = -4; a < 4; a++) {
        for (int b = -4; b < 4; b++) {
            vec2 grid = vec2(a, b);
            vec2 offset = vec2(hash(grid), hash(grid + 1.7));
            vec3 center = vec3(float(a) + 0.9 * offset.x, 0.2, float(b) + 0.9 * offset.y);

            if (length(center - vec3(4, 0.2, 0)) > 0.9) {
                float choose_mat = hash(grid + 3.5);
                Material mat;

                if (choose_mat < 0.8) {
                    vec3 albedo = vec3(hash(grid + 0.1), hash(grid + 0.2), hash(grid + 0.3)) *
                                  vec3(hash(grid + 0.4), hash(grid + 0.5), hash(grid + 0.6));
                    mat = Material(albedo, 0.0, 1.0, 0);
                } else if (choose_mat < 0.95) {
                    vec3 albedo = vec3(0.5) + 0.5 * vec3(hash(grid + 0.7), hash(grid + 0.8), hash(grid + 0.9));
                    float fuzz = 0.5 * hash(grid + 1.0);
                    mat = Material(albedo, fuzz, 1.0, 1);
                } else {
                    mat = Material(vec3(1.0), 0.0, 1.5, 2);
                }

                if (hit_sphere(center, 0.2, ro, rd, t, normal) && t < rec.t) {
                    rec.t = t;
                    rec.p = ro + t * rd;
                    rec.normal = normal;
                    rec.material = mat;
                    rec.hit = true;
                }
            }
        }
    }

    return rec;
}

vec3 ray_color(vec3 ro, vec3 rd, vec2 seed) {
    vec3 attenuation = vec3(1.0);
    for (int i = 0; i < MAX_BOUNCES; i++) {
        HitRecord rec = world_hit(ro, rd);
        if (!rec.hit) {
            vec3 unit_dir = normalize(rd);
            float t = 0.5 * (unit_dir.y + 1.0);
            vec3 c = mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            return attenuation * c;
        }

        if (rec.material.type == 0) {
            // Lambertian
            vec3 target = normalize(rec.normal + random_in_unit_sphere(seed + float(i)));
            attenuation *= rec.material.albedo;
            ro = rec.p;
            rd = target;
        }
        else if (rec.material.type == 1) {
            // Metal
            vec3 reflected = reflect(normalize(rd), rec.normal);
            rd = normalize(reflected + rec.material.fuzz * random_in_unit_sphere(seed + float(i)));
            if (dot(rd, rec.normal) <= 0.0) return vec3(0.0);
            attenuation *= rec.material.albedo;
            ro = rec.p;
        }
        else {
            // Dielectric with Schlick approx
            float refraction_ratio = dot(rd, rec.normal) < 0.0 ? 1.0 / rec.material.ref_idx : rec.material.ref_idx;
            float cos_theta = min(dot(-rd, rec.normal), 1.0);
            float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            float R0 = pow((1.0 - refraction_ratio) / (1.0 + refraction_ratio), 2.0);
            float reflect_prob = R0 + (1.0 - R0) * pow(1.0 - cos_theta, 5.0);

            if (refraction_ratio * sin_theta > 1.0 || hash(seed + float(i)) < reflect_prob) {
                rd = reflect(normalize(rd), rec.normal);
            } else {
                rd = refract(normalize(rd), rec.normal, refraction_ratio);
            }

            attenuation *= rec.material.albedo;
            ro = rec.p;
        }
    }
    return vec3(0.0);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = (fragCoord) / iResolution.xy;
    uv = 2.0 * uv - 1.0;
    uv.x *= iResolution.x / iResolution.y;

    vec3 lookfrom = vec3(13, 2, 3);
    vec3 lookat = vec3(0, 0, 0);
    vec3 vup = vec3(0, 1, 0);
    float focus_dist = 10.0;
    float vfov = radians(20.0);

    vec3 w = normalize(lookfrom - lookat);
    vec3 u = normalize(cross(vup, w));
    vec3 v = cross(w, u);

    float h = tan(vfov / 2.0);
    float viewport_height = 2.0 * h;
    float viewport_width = iResolution.x / iResolution.y * viewport_height;

    vec3 horizontal = focus_dist * viewport_width * u;
    vec3 vertical = focus_dist * viewport_height * v;
    vec3 lower_left = lookfrom - horizontal / 2.0 - vertical / 2.0 - focus_dist * w;

    vec3 col = vec3(0.0);
    int SAMPLES = 10;
    for (int i = 0; i < SAMPLES; i++) {
        vec2 jitter = (vec2(hash(fragCoord.xy + float(i)), hash(fragCoord.yx + float(i))) - 0.5) / iResolution.xy;
        vec2 uvj = uv + jitter;

        vec3 rd = normalize(lower_left + uvj.x * horizontal + uvj.y * vertical - lookfrom);
        col += ray_color(lookfrom, rd, fragCoord.xy + float(i));
    }

    col /= float(SAMPLES);
    fragColor = vec4(sqrt(col), 1.0); // gamma correction
}
