#include "ast.h"

const char *node_type_name(NodeType type) {
    switch (type) {
        case NODE_LITERAL:   return "Literal";
        case NODE_CLASS:     return "Class";
        case NODE_CONCAT:    return "Concat";
        case NODE_ALTERNATE: return "Alternate";
        case NODE_REPEAT:    return "Repeat";
        case NODE_REF:       return "Ref";
    }
    return "?";
}

Node *node_new(NodeType type, int line, int column) {
    Node *n = calloc(1, sizeof(Node));
    if (!n) { fprintf(stderr, "out of memory allocating AST node\n"); exit(1); }
    n->type = type;
    n->line = line;
    n->column = column;
    n->repeat_min = 1;
    n->repeat_max = 1;
    n->ref_index = -1;
    n->min_len = 0;
    n->max_len = 0;
    return n;
}

void node_add_child(Node *parent, Node *child) {
    if (parent->child_count + 1 > parent->child_capacity) {
        int new_cap = parent->child_capacity == 0 ? 4 : parent->child_capacity * 2;
        Node **grown = realloc(parent->children, sizeof(Node *) * (size_t)new_cap);
        if (!grown) { fprintf(stderr, "out of memory growing AST node\n"); exit(1); }
        parent->children = grown;
        parent->child_capacity = new_cap;
    }
    parent->children[parent->child_count++] = child;
}

void node_free(Node *node) {
    if (!node) return;
    for (int i = 0; i < node->child_count; i++) node_free(node->children[i]);
    free(node->children);
    free(node->masks);
    free(node);
}

static void print_indent(int indent) {
    for (int i = 0; i < indent; i++) printf("  ");
}

static void print_length_range(Node *node) {
    if (node->max_len == LEN_INF) {
        printf("  [len %d..inf]", node->min_len);
    } else {
        printf("  [len %d..%d]", node->min_len, node->max_len);
    }
}

void ast_print(Node *node, int indent, AlphabetKind kind, bool resolved) {
    if (!node) return;
    print_indent(indent);

    switch (node->type) {
        case NODE_LITERAL: {
            printf("Literal \"%.*s\"", node->text_length, node->text);
            if (resolved && node->masks) {
                printf("  ->");
                for (int i = 0; i < node->mask_count; i++) {
                    char buf[32];
                    alphabet_format_mask(kind, node->masks[i], buf, sizeof(buf));
                    printf(" %s", buf);
                }
            }
            break;
        }
        case NODE_CLASS: {
            printf("Class [%s%.*s]", node->negated ? "^" : "", node->text_length, node->text);
            if (resolved) {
                char buf[64];
                alphabet_format_mask(kind, node->class_mask, buf, sizeof(buf));
                printf("  -> %s", buf);
            }
            break;
        }
        case NODE_CONCAT:
            printf("Concat (%d)", node->child_count);
            break;
        case NODE_ALTERNATE:
            printf("Alternate (%d)", node->child_count);
            break;
        case NODE_REPEAT:
            if (node->repeat_max == REPEAT_INF) {
                printf("Repeat {%d,inf}", node->repeat_min);
            } else {
                printf("Repeat {%d,%d}", node->repeat_min, node->repeat_max);
            }
            break;
        case NODE_REF:
            printf("Ref @%.*s", node->text_length, node->text);
            if (resolved) printf("  -> motif #%d", node->ref_index);
            break;
    }

    if (resolved) print_length_range(node);
    printf("\n");

    for (int i = 0; i < node->child_count; i++) {
        ast_print(node->children[i], indent + 1, kind, resolved);
    }
}
