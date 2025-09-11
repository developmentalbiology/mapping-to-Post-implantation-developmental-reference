# GitHub Upload Guide

This guide walks you through uploading the spatial data label transfer tutorial to GitHub.

## Prerequisites

1. **GitHub Account**: Create account at [github.com](https://github.com) if you don't have one
2. **Git Installed**: Download from [git-scm.com](https://git-scm.com/) if not installed
3. **Repository Ready**: Your tutorial package should be complete and tested

## Step 1: Create GitHub Repository

1. Go to [GitHub.com](https://github.com) and sign in
2. Click the **"+"** icon in top right corner → **"New repository"**
3. Fill in repository details:
   - **Repository name**: `spatial_label_transfer_tutorial`
   - **Description**: `Advanced spatial annotation using a 5-tier assignment algorithm for spatial transcriptomics data`
   - **Visibility**: Choose Public (recommended for tutorials) or Private
   - **Initialize**: Leave "Add a README file" **unchecked** (we already have one)
   - **Add .gitignore**: Choose "R" template
   - **Choose a license**: MIT License (recommended for tutorials)
4. Click **"Create repository"**

## Step 2: Prepare Local Repository

Open terminal/command prompt and navigate to your tutorial directory:

```bash
cd "D:\xueying-work\afterPHD\Liu-lab\project1-embryo benchmarking\2024April\manuscript\NCB_rebuttal_2025_june\code\20250907_label_transfer_github_tutorial\spatial_label_transfer_tutorial"
```

## Step 3: Initialize Git Repository

```bash
# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: Spatial data label transfer tutorial

- Complete 5-tier assignment algorithm workflow
- Human CS8 embryo example with exact parameter reproduction
- Comprehensive visualization and plotting utilities
- YAML-based configuration system
- Complete documentation and tutorials"
```

## Step 4: Connect to GitHub

Replace `YOUR_USERNAME` with your actual GitHub username:

```bash
# Add GitHub repository as remote origin
git remote add origin https://github.com/YOUR_USERNAME/spatial_data_label_transfer_tutorial.git

# Set main branch
git branch -M main

# Push to GitHub
git push -u origin main
```

## Step 5: Verify Upload

1. Go to your GitHub repository: `https://github.com/YOUR_USERNAME/spatial_data_label_transfer_tutorial`
2. Check that all files are present:
   - README.md displays properly
   - All scripts in `scripts/` folder
   - Configuration files in `config/`
   - Example scripts in `examples/`
   - Documentation in `docs/`

## Step 6: Update Repository Links

After successful upload, update the repository links in README.md:

1. Replace `YOUR_USERNAME` with your actual GitHub username in:
   - `README.md` (installation instructions)
   - Any other documentation files

2. Create a new commit with the updated links:
```bash
git add README.md
git commit -m "Update repository links with actual GitHub username"
git push
```

## Step 7: Create Release (Optional but Recommended)

1. In your GitHub repository, click **"Releases"** (right sidebar)
2. Click **"Create a new release"**
3. Fill in release details:
   - **Tag version**: `v1.0.0`
   - **Release title**: `Spatial Data Label Transfer Tutorial v1.0.0`
   - **Description**: Brief summary of features and usage
4. Click **"Publish release"**

## Step 8: Add Repository Topics and Description

1. In your repository main page, click the ⚙️ (gear) icon next to "About"
2. Add **Topics**: `spatial-transcriptomics`, `rctd`, `cell-type-annotation`, `bioinformatics`, `r`, `tutorial`
3. Update **Description** if needed
4. Add **Website** URL if you have documentation hosted elsewhere

## File Structure Verification

Your uploaded repository should have this structure:

```
spatial_data_label_transfer_tutorial/
├── README.md                    # Updated with installation instructions
├── TUTORIAL.md                  
├── GITHUB_UPLOAD_GUIDE.md       # This file
├── LICENSE
├── install_r_packages.R         
├── requirements.txt             
├── config/
│   ├── default.yaml
│   └── human_CS8_spatial.yaml   # Human CS8 parameters
├── scripts/
│   ├── run_workflow.sh
│   ├── step1_apply_annotation.R
│   ├── step2_export_plots.R
│   ├── multi_tier_integrated_final.R
│   └── utils/
├── examples/
│   ├── basic_example.sh
│   └── human_embryo_example.sh  # Human CS8 example
├── docs/
│   ├── TUTORIAL.md
│   └── workflow_diagram_v2.png
├── data/
│   └── README.md                # Download instructions
└── output/                      # Will be created when users run workflow
```

## Troubleshooting

### Authentication Issues
If you get authentication errors:
```bash
# Use personal access token instead of password
# Go to GitHub Settings → Developer settings → Personal access tokens
# Generate new token with 'repo' permissions
# Use token as password when prompted
```

### File Too Large Errors
If any files are too large (>100MB):
```bash
# Remove large files and use Git LFS
git rm large_file.rds
echo "*.rds" >> .gitattributes
git lfs track "*.rds"
git add .gitattributes
git commit -m "Add Git LFS tracking for RDS files"
```

## Best Practices

1. **Keep data files external**: Don't upload large data files to GitHub
2. **Use clear commit messages**: Describe what changes you made
3. **Tag releases**: Use semantic versioning (v1.0.0, v1.1.0, etc.)
4. **Update documentation**: Keep README.md and tutorials current
5. **Test before uploading**: Ensure all scripts work correctly

## Next Steps After Upload

1. **Share the repository**: Give the GitHub URL to users
2. **Monitor issues**: Users may report bugs or ask questions
3. **Accept contributions**: Consider pull requests for improvements
4. **Create documentation website**: Consider GitHub Pages for better documentation

Your tutorial is now publicly available and ready for users to download and use!